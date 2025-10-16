from __future__ import annotations

import asyncio
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from backend.llm.annotator import Annotator
from config.settings import get_settings
from gpt_cell_annotator.scanpy import (
    MARKER_DB_COLUMNS,
    BatchOptions,
    DiskAnnotationCache,
    GuardrailConfig,
    annotate_anndata,
    annotate_anndata_async,
    annotate_from_markers,
    annotate_rank_genes,
    main as scanpy_main,
    validate_anndata,
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
    matrix = np.zeros((6, 4), dtype=float)
    matrix[:3, 0] = 5.0  # cluster 0 expresses gene0 (MS4A1)
    matrix[3:, 1] = 5.0  # cluster 1 expresses gene1 (CD3E)
    obs = pd.DataFrame(
        {"cluster": ["0", "0", "0", "1", "1", "1"]},
        index=[f"cell_{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=["MS4A1", "CD3E", "GNLY", "LYZ"])
    adata = ad.AnnData(matrix, obs=obs, var=var)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def test_annotate_anndata_adds_columns(monkeypatch):
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    adata = _build_adata()
    marker_db = _build_mock_marker_db()
    annotator = Annotator(settings=settings)

    result = annotate_anndata(
        adata,
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
        annotator=annotator,
    )

    assert "gptca_label" in result.adata.obs.columns
    assert "gptca_proposed_label" in result.adata.obs.columns
    labels = {label for label in result.adata.obs["gptca_label"] if label}
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
    assert result.report.summary.total_clusters == 2
    assert len(result.report.clusters) == 2
    assert result.stats["total_clusters"] == 2


def test_cli_roundtrip(tmp_path: Path, monkeypatch):
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    adata = _build_adata()
    input_path = tmp_path / "demo.h5ad"
    adata.write(input_path)

    marker_db = _build_mock_marker_db()

    monkeypatch.setattr(
        "gpt_cell_annotator.scanpy._load_marker_db",
        lambda marker_db_arg, marker_db_path, cache=True: marker_db,
    )

    summary_csv = tmp_path / "summary.csv"
    summary_json = tmp_path / "summary.json"
    stats_json = tmp_path / "stats.json"
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
            "--stats-json",
            str(stats_json),
            "--offline",
        ]
    )

    assert exit_code == 0
    annotated = ad.read_h5ad(output_path)
    assert "gptca_label" in annotated.obs.columns
    assert "gptca_proposed_label" in annotated.obs.columns
    assert summary_csv.exists()
    assert summary_json.exists()
    assert stats_json.exists()


def test_annotate_anndata_async_matches_sync():
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    marker_db = _build_mock_marker_db()

    sync_result = annotate_anndata(
        _build_adata(),
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
    )

    async_result = asyncio.run(
        annotate_anndata_async(
            _build_adata(),
            "cluster",
            species="Homo sapiens",
            marker_db=marker_db,
        )
    )

    assert async_result.report.summary.total_clusters == sync_result.report.summary.total_clusters
    assert async_result.stats["total_clusters"] == sync_result.stats["total_clusters"]
    assert "gptca_label" in async_result.adata.obs.columns


def test_annotate_from_markers(tmp_path: Path):
    marker_db = _build_mock_marker_db()
    markers = {"0": ["MS4A1", "CD79A"], "1": ["CD3E", "CD3D"]}
    result = annotate_from_markers(
        markers,
        species="Homo sapiens",
        marker_db=marker_db,
        batch_options=BatchOptions(size=1, concurrency=1),
    )

    assert result.report.summary.total_clusters == 2
    labels = {cluster.annotation["primary_label"] for cluster in result.report.clusters}
    assert "Unknown or Novel" not in labels


def test_annotate_rank_genes_wrapper():
    marker_db = _build_mock_marker_db()
    rankings = {
        "0": ["MS4A1", "CD79A", "CD74"],
        "1": ["CD3E", "CD3D", "CD2"],
    }
    result = annotate_rank_genes(
        rankings,
        species="Homo sapiens",
        marker_db=marker_db,
    )
    assert result.report.summary.total_clusters == 2


def test_validate_anndata_guardrail_override():
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    adata = _build_adata()
    adata.obs["gpt_label"] = ["B cell"] * 6
    marker_db = _build_mock_marker_db()

    report_default = validate_anndata(
        adata,
        "cluster",
        species="Homo sapiens",
        label_column="gpt_label",
        marker_db=marker_db,
    )
    assert report_default.summary.flagged_clusters <= 2

    report_strict = validate_anndata(
        adata,
        "cluster",
        species="Homo sapiens",
        label_column="gpt_label",
        marker_db=marker_db,
        guardrails=GuardrailConfig(min_marker_overlap=5),
    )
    assert report_strict.summary.flagged_clusters >= report_default.summary.flagged_clusters


def test_disk_cache_reduces_batches(tmp_path: Path, monkeypatch):
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    marker_db = _build_mock_marker_db()

    call_counter = {"count": 0}

    cache_dir = tmp_path / "cache"
    cache = DiskAnnotationCache(cache_dir)

    monkeypatch.setattr(
        "gpt_cell_annotator.scanpy.Annotator.annotate_batch",
        lambda self, payload, context: call_counter.__setitem__("count", call_counter["count"] + 1)
        or {str(item["cluster_id"]): {"primary_label": "B cell"} for item in payload},
    )

    result_first = annotate_anndata(
        _build_adata(),
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
        annotation_cache=cache,
        batch_options=BatchOptions(size=1, concurrency=1),
    )
    assert result_first.stats["llm_batches"] == 2

    call_counter["count"] = 0
    result_second = annotate_anndata(
        _build_adata(),
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
        annotation_cache=cache,
        batch_options=BatchOptions(size=1, concurrency=1),
    )
    assert call_counter["count"] == 0
    assert result_second.stats["cache_hits"] == 2
