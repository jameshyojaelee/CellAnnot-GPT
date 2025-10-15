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
- evidence_score: heuristic label capturing strength of the evidence text.

The data can be materialized as Parquet (columnar analytics) and/or SQLite for
lightweight querying inside the application.
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import sqlite3
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import httpx
import pandas as pd
import yaml

logger = logging.getLogger("cellannot.data_ingest")

NORMALIZED_COLUMNS = [
    "source",
    "cell_type",
    "ontology_id",
    "gene_symbol",
    "species",
    "tissue",
    "evidence",
    "reference",
    "evidence_score",
]

EVIDENCE_SCORE_COLUMN = "evidence_score"


@dataclass
class SourceConfig:
    """Configuration for a marker gene source."""

    name: str
    fmt: str
    parser: Callable[[bytes, str], pd.DataFrame]
    metadata: dict[str, str] = field(default_factory=dict)
    url: str | None = None
    local_path: Path | None = None
    required_columns: Sequence[str] | None = None
    checksum: str | None = None
    checksum_algorithm: str = "sha256"


class SourceResolutionError(RuntimeError):
    """Raised when a source cannot be resolved from local or remote locations."""


class ChecksumMismatchError(RuntimeError):
    """Raised when the downloaded payload does not match the expected checksum."""


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
    normalized[EVIDENCE_SCORE_COLUMN] = ""
    return normalized.reindex(columns=NORMALIZED_COLUMNS, fill_value="")


def parse_cellmarker(payload: bytes, source: str) -> pd.DataFrame:
    """Parse CellMarker CSV export into the normalized schema."""

    df = pd.read_csv(io.BytesIO(payload))
    expected_cols = {"cell_type", "gene_symbol", "tissue_type", "species", "pubmed_id"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"{source} is missing columns: {', '.join(sorted(missing))}")

    def _pubmed_url(pubmed_id: Any) -> str:
        if pd.notna(pubmed_id):
            return f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}"
        return ""

    references = df["pubmed_id"].map(_pubmed_url)

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
    normalized[EVIDENCE_SCORE_COLUMN] = ""
    return normalized.reindex(columns=NORMALIZED_COLUMNS, fill_value="")


def parse_curated_json(payload: bytes, source: str) -> pd.DataFrame:
    """Parse curated JSON snippets into the normalized schema."""

    records = json.loads(payload.decode("utf-8"))
    if not isinstance(records, list):
        raise ValueError(f"{source} must be a list of marker objects")

    normalized: list[dict[str, str]] = []
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


def compute_checksum(payload: bytes, algorithm: str) -> str:
    """Return the hexadecimal digest for the given payload."""

    try:
        hasher = getattr(hashlib, algorithm)
    except AttributeError as exc:
        raise ValueError(f"Unsupported checksum algorithm '{algorithm}'") from exc
    return hasher(payload).hexdigest()


def classify_evidence_strength(evidence: str) -> str:
    """Heuristically score evidence strings as low/medium/high."""

    if not evidence:
        return "low"

    text = evidence.lower()
    score = 0
    if any(token in text for token in ("logfc", "fold change", "lfc", "log2fc", "effect size")):
        score += 2
    if any(token in text for token in ("p-value", "p<", "adjusted p", "padj", "fdr", "q-value")):
        score += 2
    if any(token in text for token in ("single-cell", "scrna", "rna-seq", "bulk", "proteomics")):
        score += 1
    if any(token in text for token in ("validated", "marker", "signature", "score")):
        score += 1

    if score >= 4:
        return "high"
    if score >= 2:
        return "medium"
    return "low"


def ensure_evidence_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Populate the evidence_score column if missing."""

    if EVIDENCE_SCORE_COLUMN not in df.columns:
        df[EVIDENCE_SCORE_COLUMN] = df["evidence"].fillna("").map(classify_evidence_strength)
    else:
        df[EVIDENCE_SCORE_COLUMN] = (
            df[EVIDENCE_SCORE_COLUMN]
            .fillna("")
            .where(
                df[EVIDENCE_SCORE_COLUMN].astype(bool),
                df["evidence"].fillna("").map(classify_evidence_strength),
            )
        )
    return df


PARSER_REGISTRY: dict[str, Callable[[bytes, str], pd.DataFrame]] = {
    "panglaodb": parse_panglaodb,
    "cellmarker": parse_cellmarker,
    "curated_json": parse_curated_json,
    "csv": parse_panglaodb,  # default CSV structure
    "json": parse_curated_json,
}


class MarkerDataLoader:
    """Download, normalize, and persist marker gene knowledge sources."""

    def __init__(
        self,
        sources: Iterable[SourceConfig],
        storage_dir: Path,
        *,
        parquet_path: Path | None = None,
        sqlite_path: Path | None = None,
        http_client: httpx.Client | None = None,
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

    def _fetch(
        self,
        source: SourceConfig,
        *,
        local_only: bool = False,
        enforce_checksums: bool = False,
    ) -> bytes:
        """Return raw bytes for a source, preferring local assets when present."""

        resolution_path = None
        if source.local_path and source.local_path.exists():
            resolution_path = str(source.local_path)
            payload = source.local_path.read_bytes()
        else:
            if local_only:
                raise SourceResolutionError(
                    f"Local file for source '{source.name}' not found and local-only mode enabled."
                )
            if not source.url:
                raise SourceResolutionError(f"Source '{source.name}' missing URL and local path.")
            try:
                response = self._http_client.get(
                    str(source.url),
                    follow_redirects=True,
                )
                response.raise_for_status()
            except httpx.HTTPError as exc:
                msg = f"Failed to download '{source.name}' from {source.url}: {exc}"
                raise SourceResolutionError(msg) from exc
            payload = response.content
            resolution_path = source.url

        if source.checksum:
            digest = compute_checksum(payload, source.checksum_algorithm)
            if enforce_checksums and digest.lower() != source.checksum.lower():
                raise ChecksumMismatchError(
                    f"{source.name} checksum mismatch: expected {source.checksum}, got {digest}"
                )
            if digest.lower() != (source.checksum or "").lower():
                logger.debug("Checksum for %s (computed=%s)", source.name, digest)
        logger.debug("Resolved source %s via %s", source.name, resolution_path)
        return payload

    def _validate_schema(self, df: pd.DataFrame, source: SourceConfig) -> None:
        """Ensure the normalized dataset contains required columns."""

        required: set[str]
        if source.required_columns:
            required = {col for col in source.required_columns}
        else:
            required = set(NORMALIZED_COLUMNS) - {EVIDENCE_SCORE_COLUMN}

        missing = required - set(df.columns)
        if missing:
            missing_columns = ", ".join(sorted(missing))
            raise ValueError(
                f"{source.name} missing required columns after parsing: {missing_columns}"
            )

    def _normalise_frame(self, df: pd.DataFrame, source: SourceConfig) -> pd.DataFrame:
        """Apply common cleanup to parsed frames."""

        self._validate_schema(df, source)
        df = ensure_evidence_scores(df)
        if source.metadata.get("version"):
            df["source_version"] = source.metadata["version"]
        else:
            df["source_version"] = ""
        return df

    def load_all(
        self,
        *,
        local_only: bool = False,
        enforce_checksums: bool = False,
    ) -> pd.DataFrame:
        frames: list[pd.DataFrame] = []
        resolved_sources = 0
        resolution_errors: list[str] = []

        for cfg in self.sources:
            try:
                payload = self._fetch(
                    cfg,
                    local_only=local_only,
                    enforce_checksums=enforce_checksums,
                )
            except ChecksumMismatchError as exc:
                logger.error("Checksum verification failed for %s: %s", cfg.name, exc)
                raise
            except SourceResolutionError as exc:
                resolution_errors.append(f"{cfg.name}: {exc}")
                logger.warning("%s", exc)
                continue

            resolved_sources += 1
            try:
                df = cfg.parser(payload, cfg.name)
            except Exception as exc:  # - bubble up parser issues with context
                raise ValueError(f"{cfg.name} parser failed: {exc}") from exc

            df = self._normalise_frame(df, cfg)
            frames.append(df)

        if resolved_sources == 0:
            joined = "; ".join(resolution_errors) or "no sources configured"
            raise SourceResolutionError(f"Failed to resolve any marker sources ({joined})")

        if not frames:
            raise ValueError("Resolved sources but no marker records were produced")

        combined = pd.concat(frames, ignore_index=True)
        combined = combined.drop_duplicates(
            subset=["source", "cell_type", "gene_symbol", "species"]
        )
        # ensure column order for downstream consumers
        ordered_cols = [col for col in NORMALIZED_COLUMNS if col in combined.columns]
        remaining_cols = [col for col in combined.columns if col not in ordered_cols]
        combined = combined.loc[:, ordered_cols + remaining_cols]
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

    def run(
        self,
        *,
        write_parquet: bool = True,
        write_sqlite: bool = True,
        local_only: bool = False,
        enforce_checksums: bool = False,
    ) -> pd.DataFrame:
        try:
            df = self.load_all(local_only=local_only, enforce_checksums=enforce_checksums)
            if write_parquet:
                self.write_to_parquet(df)
            if write_sqlite:
                self.write_to_sqlite(df)
            return df
        finally:
            self.close()


def load_sources_from_yaml(path: Path) -> list[SourceConfig]:
    if not path.exists():
        raise FileNotFoundError(f"Marker source config not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    sources_cfg = data.get("sources", [])
    sources: list[SourceConfig] = []
    for entry in sources_cfg:
        name = entry["name"]
        fmt = entry.get("fmt", "csv").lower()
        parser_key = entry.get("parser", fmt)
        parser = PARSER_REGISTRY.get(parser_key)
        if parser is None:
            raise ValueError(f"Unknown parser '{parser_key}' for source '{name}'")
        local_path = entry.get("local_path")
        source = SourceConfig(
            name=name,
            fmt=fmt,
            parser=parser,
            metadata=entry.get("metadata", {}),
            url=entry.get("url"),
            local_path=Path(local_path) if local_path else None,
            required_columns=entry.get("required_columns"),
            checksum=entry.get("checksum"),
            checksum_algorithm=entry.get("checksum_algorithm", "sha256"),
        )
        sources.append(source)
    return sources


def default_sources(config_path: Path | None = None) -> list[SourceConfig]:
    if config_path:
        return load_sources_from_yaml(config_path)
    return []


__all__ = [
    "MarkerDataLoader",
    "SourceConfig",
    "parse_panglaodb",
    "parse_cellmarker",
    "parse_curated_json",
    "load_sources_from_yaml",
    "default_sources",
    "NORMALIZED_COLUMNS",
]
