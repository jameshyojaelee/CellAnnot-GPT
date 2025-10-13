from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from backend.data_ingest.marker_loader import (
    MarkerDataLoader,
    SourceConfig,
    NORMALIZED_COLUMNS,
    parse_cellmarker,
    parse_curated_json,
    parse_panglaodb,
)


class FakeResponse:
    def __init__(self, content: bytes) -> None:
        self.content = content
        self.status_code = 200

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise ValueError("HTTP error")


class FakeClient:
    def __init__(self, payloads: dict[str, bytes]) -> None:
        self._payloads = payloads

    def get(self, url: str) -> FakeResponse:
        try:
            return FakeResponse(self._payloads[url])
        except KeyError as exc:
            raise ValueError(f"Unexpected URL requested: {url}") from exc


@pytest.fixture()
def sample_payloads() -> dict[str, bytes]:
    panglao_csv = (
        "cell_type,gene,organ,species,evidence,reference\n"
        "CD4 T cell,CCR7,Blood,Homo sapiens,logFC=1.2,https://panglaodb.se\n"
    )
    cellmarker_csv = (
        "cell_type,gene_symbol,tissue_type,species,pubmed_id\n"
        "B cell,MS4A1,Blood,Homo sapiens,12345678\n"
    )
    curated_json = json.dumps(
        [
            {
                "cell_type": "Alveolar macrophage",
                "ontology_id": "CL:0000583",
                "gene_symbol": "MARCO",
                "species": "Homo sapiens",
                "tissue": "Lung",
                "evidence": "Expressed in scRNA-seq study",
                "reference": "https://doi.org/10.1038/example",
            }
        ]
    )
    return {
        "https://storage.example.org/panglaodb_markers.csv": panglao_csv.encode(),
        "https://storage.example.org/cellmarker_markers.csv": cellmarker_csv.encode(),
        "https://storage.example.org/curated_markers.json": curated_json.encode(),
    }


def test_marker_loader_normalizes_sources(tmp_path: Path, sample_payloads: dict[str, bytes]) -> None:
    client = FakeClient(sample_payloads)
    sources = [
        SourceConfig(
            name="PanglaoDB",
            url="https://storage.example.org/panglaodb_markers.csv",
            fmt="csv",
            parser=parse_panglaodb,
        ),
        SourceConfig(
            name="CellMarker",
            url="https://storage.example.org/cellmarker_markers.csv",
            fmt="csv",
            parser=parse_cellmarker,
        ),
        SourceConfig(
            name="CuratedLiterature",
            url="https://storage.example.org/curated_markers.json",
            fmt="json",
            parser=parse_curated_json,
        ),
    ]

    loader = MarkerDataLoader(
        sources,
        storage_dir=tmp_path,
        parquet_path=tmp_path / "markers.parquet",
        sqlite_path=tmp_path / "markers.sqlite",
        http_client=client,
    )

    df = loader.run(write_parquet=True, write_sqlite=True)
    assert set(df.columns) == set(NORMALIZED_COLUMNS)
    assert len(df) == 3

    parquet_path = tmp_path / "markers.parquet"
    sqlite_path = tmp_path / "markers.sqlite"
    assert parquet_path.exists()
    assert sqlite_path.exists()

    roundtrip_df = pd.read_parquet(parquet_path)
    assert len(roundtrip_df) == 3

    # Validate SQLite contents
    import sqlite3

    with sqlite3.connect(sqlite_path) as conn:
        count = conn.execute("SELECT COUNT(*) FROM cell_markers").fetchone()[0]
    assert count == 3


def test_parser_error_when_columns_missing() -> None:
    # Panglao parser expects a column called gene
    bad_payload = "cell_type,organ,species\nfoo,bar,baz\n".encode()
    with pytest.raises(ValueError):
        parse_panglaodb(bad_payload, "PanglaoDB")


def test_curated_parser_requires_list() -> None:
    bad_payload = json.dumps({"cell_type": "T cell"}).encode()
    with pytest.raises(ValueError):
        parse_curated_json(bad_payload, "CuratedLit")
