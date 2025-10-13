"""Utilities for building the marker-gene knowledge base.

The ingestion pipeline normalizes multiple public and curated sources into a
unified tabular schema with the following columns:

- source: canonical name of the upstream dataset.
- cell_type: human-readable cell type label.
- ontology_id: optional identifier from Cell Ontology or similar.
- gene_symbol: HGNC-aligned gene symbol.
- species: species in which the marker is reported (e.g. "Homo sapiens").
- tissue: tissue or compartment context when available.
- evidence: free-text evidence string (score, logFC, citation snippet, etc.).
- reference: URL or citation to the supporting resource.

The data can be materialized as Parquet (columnar analytics) and/or SQLite for
lightweight querying inside the application.
"""

from __future__ import annotations

import io
import json
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

import httpx
import pandas as pd

NORMALIZED_COLUMNS = [
    "source",
    "cell_type",
    "ontology_id",
    "gene_symbol",
    "species",
    "tissue",
    "evidence",
    "reference",
]


@dataclass
class SourceConfig:
    """Configuration for a marker gene source."""

    name: str
    url: str
    fmt: str
    parser: Callable[[bytes, str], pd.DataFrame]
    metadata: Dict[str, str] = field(default_factory=dict)
    local_path: Optional[Path] = None


def parse_panglaodb(payload: bytes, source: str) -> pd.DataFrame:
    """Parse PanglaoDB CSV export into the normalized schema."""

    df = pd.read_csv(io.BytesIO(payload))
    expected_cols = {"cell_type", "gene", "organ", "species", "evidence", "reference"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"{source} is missing columns: {', '.join(sorted(missing))}")

    normalized = pd.DataFrame(
        {
            "source": source,
            "cell_type": df["cell_type"],
            "ontology_id": df.get("ontology_id", ""),
            "gene_symbol": df["gene"],
            "species": df["species"],
            "tissue": df["organ"],
            "evidence": df["evidence"],
            "reference": df["reference"],
        }
    )
    return normalized[NORMALIZED_COLUMNS]


def parse_cellmarker(payload: bytes, source: str) -> pd.DataFrame:
    """Parse CellMarker CSV export into the normalized schema."""

    df = pd.read_csv(io.BytesIO(payload))
    expected_cols = {"cell_type", "gene_symbol", "tissue_type", "species", "pubmed_id"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"{source} is missing columns: {', '.join(sorted(missing))}")

    references = df["pubmed_id"].map(lambda x: f"https://pubmed.ncbi.nlm.nih.gov/{x}" if pd.notna(x) else "")

    normalized = pd.DataFrame(
        {
            "source": source,
            "cell_type": df["cell_type"],
            "ontology_id": df.get("ontology_id", ""),
            "gene_symbol": df["gene_symbol"],
            "species": df["species"],
            "tissue": df["tissue_type"],
            "evidence": df.get("evidence", ""),
            "reference": references,
        }
    )
    return normalized[NORMALIZED_COLUMNS]


def parse_curated_json(payload: bytes, source: str) -> pd.DataFrame:
    """Parse curated JSON snippets into the normalized schema."""

    records = json.loads(payload.decode("utf-8"))
    if not isinstance(records, list):
        raise ValueError(f"{source} must be a list of marker objects")

    normalized: List[Dict[str, str]] = []
    for record in records:
        normalized.append(
            {
                "source": source,
                "cell_type": record.get("cell_type", ""),
                "ontology_id": record.get("ontology_id", ""),
                "gene_symbol": record.get("gene_symbol", ""),
                "species": record.get("species", ""),
                "tissue": record.get("tissue", ""),
                "evidence": record.get("evidence", ""),
                "reference": record.get("reference", ""),
            }
        )

    df = pd.DataFrame(normalized, columns=NORMALIZED_COLUMNS)
    return df


class MarkerDataLoader:
    """Download, normalize, and persist marker gene knowledge sources."""

    def __init__(
        self,
        sources: Iterable[SourceConfig],
        storage_dir: Path,
        *,
        parquet_path: Optional[Path] = None,
        sqlite_path: Optional[Path] = None,
        http_client: Optional[httpx.Client] = None,
    ) -> None:
        self.sources = list(sources)
        self.storage_dir = storage_dir
        self.storage_dir.mkdir(parents=True, exist_ok=True)
        self.parquet_path = parquet_path or (self.storage_dir / "marker_db.parquet")
        self.sqlite_path = sqlite_path or (self.storage_dir / "marker_db.sqlite")
        self._http_client = http_client or httpx.Client(timeout=30.0)
        self._owns_client = http_client is None

    def close(self) -> None:
        if self._owns_client:
            self._http_client.close()

    def _fetch(self, source: SourceConfig) -> bytes:
        if source.local_path and source.local_path.exists():
            return source.local_path.read_bytes()
        response = self._http_client.get(source.url)
        response.raise_for_status()
        return response.content

    def load_all(self) -> pd.DataFrame:
        frames: List[pd.DataFrame] = []
        for cfg in self.sources:
            payload = self._fetch(cfg)
            df = cfg.parser(payload, cfg.name)
            missing = set(NORMALIZED_COLUMNS) - set(df.columns)
            if missing:
                raise ValueError(f"{cfg.name} parser returned missing columns: {missing}")
            frames.append(df)
        if not frames:
            raise ValueError("No data sources configured")
        combined = pd.concat(frames, ignore_index=True)
        combined = combined.drop_duplicates(subset=["source", "cell_type", "gene_symbol", "species"])
        return combined

    def write_to_parquet(self, df: pd.DataFrame) -> Path:
        self.parquet_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(self.parquet_path, index=False)
        return self.parquet_path

    def write_to_sqlite(self, df: pd.DataFrame, table_name: str = "cell_markers") -> Path:
        self.sqlite_path.parent.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(self.sqlite_path) as conn:
            df.to_sql(table_name, conn, if_exists="replace", index=False)
        return self.sqlite_path

    def run(self, *, write_parquet: bool = True, write_sqlite: bool = True) -> pd.DataFrame:
        try:
            df = self.load_all()
            if write_parquet:
                self.write_to_parquet(df)
            if write_sqlite:
                self.write_to_sqlite(df)
            return df
        finally:
            self.close()


def default_sources() -> List[SourceConfig]:
    """Return SourceConfig entries for default public datasets."""

    return [
        SourceConfig(
            name="PanglaoDB",
            url="https://storage.example.org/panglaodb_markers.csv",
            fmt="csv",
            parser=parse_panglaodb,
            metadata={"license": "CC0"},
        ),
        SourceConfig(
            name="CellMarker",
            url="https://storage.example.org/cellmarker_markers.csv",
            fmt="csv",
            parser=parse_cellmarker,
            metadata={"license": "CC BY 4.0"},
        ),
        SourceConfig(
            name="CuratedLiterature",
            url="https://storage.example.org/curated_markers.json",
            fmt="json",
            parser=parse_curated_json,
            metadata={"description": "Hand-curated niche markers"},
        ),
    ]


__all__ = [
    "MarkerDataLoader",
    "SourceConfig",
    "default_sources",
    "parse_panglaodb",
    "parse_cellmarker",
    "parse_curated_json",
    "NORMALIZED_COLUMNS",
]
