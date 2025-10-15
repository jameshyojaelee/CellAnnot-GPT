from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


class _DummyCompletions:
    def __init__(self):
        self.create = lambda *args, **kwargs: None


class _DummyOpenAI:
    def __init__(self, *args, **kwargs):
        self.chat = SimpleNamespace(completions=_DummyCompletions())


sys.modules.setdefault("openai", SimpleNamespace(OpenAI=_DummyOpenAI))

from backend.llm.annotator import Annotator
from config.settings import Settings, get_settings
from gpt_cell_annotator.scanpy import (
    MARKER_DB_COLUMNS,
    annotate_anndata,
    main as scanpy_main,
)


def _build_mock_marker_db() -> pd.DataFrame:
    records = [
        {
            "source": "Demo",
            "cell_type": "B cell",
            "ontology_id": "CL:0000236",
            "gene_symbol": "MS4A1",
            "species": "Homo sapiens",
            "tissue": "Blood",
            "evidence": "Demo",
            "reference": "",
            "evidence_score": "high",
        },
        {
            "source": "Demo",
            "cell_type": "CD4 T cell",
            "ontology_id": "CL:0000624",
            "gene_symbol": "CD3E",
            "species": "Homo sapiens",
            "tissue": "Blood",
            "evidence": "Demo",
            "reference": "",
            "evidence_score": "high",
        },
    ]
    return pd.DataFrame.from_records(records, columns=MARKER_DB_COLUMNS)


def _build_adata() -> ad.AnnData:
    X = np.zeros((6, 4), dtype=float)
    X[:3, 0] = 5.0  # cluster 0 expresses gene0 (MS4A1)
    X[3:, 1] = 5.0  # cluster 1 expresses gene1 (CD3E)
    obs = pd.DataFrame(
        {"cluster": ["0", "0", "0", "1", "1", "1"]},
        index=[f"cell_{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=["MS4A1", "CD3E", "GNLY", "LYZ"])
    adata = ad.AnnData(X, obs=obs, var=var)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def test_annotate_anndata_adds_columns(monkeypatch):
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    adata = _build_adata()
    marker_db = _build_mock_marker_db()
    settings = Settings(openai_api_key="", openai_requests_per_minute=0)
    annotator = Annotator(settings=settings)

    updated, report = annotate_anndata(
        adata,
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
        annotator=annotator,
    )

    assert "gptca_label" in updated.obs.columns
    assert "gptca_proposed_label" in updated.obs.columns
    labels = {label for label in updated.obs["gptca_label"] if label}
    assert labels  # at least one label populated
    allowed = {
        "B cell",
        "CD4 T cell",
        "CD8 T cell",
        "NK cell",
        "Monocyte",
        "Platelet",
        "Erythrocyte",
        "Unknown or Novel",
    }
    assert labels <= allowed
    assert report["summary"]["total_clusters"] == 2
    assert len(report["clusters"]) == 2


def test_cli_roundtrip(tmp_path: Path, monkeypatch):
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    adata = _build_adata()
    input_path = tmp_path / "demo.h5ad"
    adata.write(input_path)

    marker_db = _build_mock_marker_db()

    monkeypatch.setattr(
        "gpt_cell_annotator.scanpy._load_marker_db",
        lambda marker_db_arg, marker_db_path: marker_db,
    )

    summary_csv = tmp_path / "summary.csv"
    summary_json = tmp_path / "summary.json"
    output_path = tmp_path / "annotated.h5ad"

    exit_code = scanpy_main(
        [
            "annotate",
            str(input_path),
            "--cluster-key",
            "cluster",
            "--species",
            "Homo sapiens",
            "--output",
            str(output_path),
            "--summary-csv",
            str(summary_csv),
            "--summary-json",
            str(summary_json),
        ]
    )

    assert exit_code == 0
    annotated = ad.read_h5ad(output_path)
    assert "gptca_label" in annotated.obs.columns
    assert "gptca_proposed_label" in annotated.obs.columns
    assert summary_csv.exists()
    assert summary_json.exists()
