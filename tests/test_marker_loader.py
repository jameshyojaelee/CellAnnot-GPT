from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from backend.data_ingest.marker_loader import (
    NORMALIZED_COLUMNS,
    ChecksumMismatchError,
    MarkerDataLoader,
    SourceConfig,
    SourceResolutionError,
    load_sources_from_yaml,
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

    def get(self, url: str, follow_redirects: bool = False) -> FakeResponse:
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


def test_marker_loader_normalizes_sources(
    tmp_path: Path, sample_payloads: dict[str, bytes]
) -> None:
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
    assert set(NORMALIZED_COLUMNS).issubset(df.columns)
    assert len(df) == 3
    assert "source_version" in df.columns
    scores = dict(
        zip(
            df["cell_type"],
            df["evidence_score"],
            strict=False,
        )
    )
    assert scores["CD4 T cell"] == "medium"

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
    bad_payload = b"cell_type,organ,species\nfoo,bar,baz\n"
    with pytest.raises(ValueError):
        parse_panglaodb(bad_payload, "PanglaoDB")


def test_curated_parser_requires_list() -> None:
    bad_payload = json.dumps({"cell_type": "T cell"}).encode()
    with pytest.raises(ValueError):
        parse_curated_json(bad_payload, "CuratedLit")


def test_load_sources_from_yaml_local_only(tmp_path: Path) -> None:
    panglao_path = tmp_path / "panglao.csv"
    panglao_path.write_text(
        "cell_type,gene,organ,species,evidence,reference\n"
        "B cell,MS4A1,Blood,Homo sapiens,known,ref\n",
        encoding="utf-8",
    )
    cellmarker_path = tmp_path / "cellmarker.csv"
    cellmarker_path.write_text(
        "cell_type,gene_symbol,tissue_type,species,pubmed_id\n"
        "T cell,CD3E,Blood,Homo sapiens,123456\n",
        encoding="utf-8",
    )
    config_path = tmp_path / "sources.yaml"
    config_text = (
        "sources:\n"
        f"  - name: PanglaoLocal\n    fmt: csv\n    parser: panglaodb\n"
        f"    local_path: {panglao_path}\n"
        f"  - name: CellMarkerLocal\n    fmt: csv\n    parser: cellmarker\n"
        f"    local_path: {cellmarker_path}\n"
    )
    config_path.write_text(config_text, encoding="utf-8")

    sources = load_sources_from_yaml(config_path)
    loader = MarkerDataLoader(sources, storage_dir=tmp_path)
    df = loader.run(local_only=True)
    assert len(df) == 2
    assert set(NORMALIZED_COLUMNS).issubset(df.columns)
    assert df["evidence_score"].isin({"low", "medium", "high"}).all()


def test_schema_validation_catches_missing_columns(tmp_path: Path) -> None:
    def broken_parser(payload: bytes, source: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "source": [source],
                "cell_type": ["foo"],
                "species": ["Homo sapiens"],
            }
        )

    sources = [
        SourceConfig(
            name="Broken",
            fmt="csv",
            parser=broken_parser,
            url="https://example.org/broken.csv",
        )
    ]
    client = FakeClient({"https://example.org/broken.csv": b"dummy"})
    loader = MarkerDataLoader(sources, storage_dir=tmp_path, http_client=client)

    with pytest.raises(ValueError, match="missing required columns"):
        loader.run()


def test_checksum_enforcement(tmp_path: Path) -> None:
    content = (
        b"cell_type,gene,organ,species,evidence,reference\n"
        b"B cell,MS4A1,Blood,Homo sapiens,Validated marker,https://example.org\n"
    )
    sources = [
        SourceConfig(
            name="PanglaoChecksum",
            fmt="csv",
            parser=parse_panglaodb,
            url="https://example.org/panglao.csv",
            checksum="deadbeef",
        )
    ]
    client = FakeClient({"https://example.org/panglao.csv": content})
    loader = MarkerDataLoader(sources, storage_dir=tmp_path, http_client=client)

    with pytest.raises(ChecksumMismatchError):
        loader.run(enforce_checksums=True)


def test_source_resolution_error(tmp_path: Path) -> None:
    sources = [
        SourceConfig(
            name="Missing",
            fmt="csv",
            parser=parse_panglaodb,
            local_path=tmp_path / "does_not_exist.csv",
        )
    ]
    loader = MarkerDataLoader(sources, storage_dir=tmp_path)

    with pytest.raises(SourceResolutionError):
        loader.run(local_only=True)


def test_fixture_sources_yaml_round_trip(tmp_path: Path) -> None:
    fixtures_yaml = Path(__file__).parent / "fixtures" / "marker_sources" / "sources.yaml"
    sources = load_sources_from_yaml(fixtures_yaml)
    loader = MarkerDataLoader(sources, storage_dir=tmp_path)
    df = loader.run(local_only=True)
    assert len(df) == 4
    assert df["evidence_score"].isin({"low", "medium", "high"}).all()
    assert set(df["source"]).issuperset({"PanglaoFixture", "CellMarkerFixture"})
