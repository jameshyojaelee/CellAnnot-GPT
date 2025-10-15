from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

import pandas as pd

class _DummyOpenAI:
    def __init__(self, *_, **__):
        self.chat = SimpleNamespace(completions=SimpleNamespace(create=lambda **kwargs: None))


sys.modules.setdefault("openai", SimpleNamespace(OpenAI=_DummyOpenAI))


from evaluation.benchmark_runner import (  # noqa: E402
    BenchmarkResult,
    load_and_run,
    run_gpt_benchmark,
    run_marker_overlap_baseline,
)
from evaluation.report_templates import render_dataset_report, render_sparkline_csv
from config.settings import get_settings

MARKER_DB = pd.DataFrame(
    [
        {
            "source": "Demo",
            "cell_type": "B cell",
            "ontology_id": "CL:0000236",
            "gene_symbol": "MS4A1",
            "species": "Homo sapiens",
            "tissue": "Blood",
            "evidence": "",
            "reference": "",
            "evidence_score": "high",
        },
        {
            "source": "Demo",
            "cell_type": "T cell",
            "ontology_id": "CL:0000084",
            "gene_symbol": "CD3D",
            "species": "Homo sapiens",
            "tissue": "Blood",
            "evidence": "",
            "reference": "",
            "evidence_score": "high",
        },
    ]
)

class DummyAnnotator:
    def annotate_cluster(self, payload, dataset_context=None):  # type: ignore[override]
        cluster_id = payload["cluster_id"]
        if cluster_id == "0":
            return {
                "primary_label": "B cell",
                "ontology_id": "CL:0000236",
                "alternatives": [{"label": "Plasma cell", "reason": "Overlap"}],
                "rationale": "MS4A1 detected",
            }
        return {
            "primary_label": "T cell",
            "ontology_id": "CL:0000084",
            "alternatives": [],
            "rationale": "Default",
        }

    def annotate_batch(self, clusters, dataset_context=None):  # type: ignore[override]
        return {
            str(cluster["cluster_id"]): self.annotate_cluster(cluster, dataset_context)
            for cluster in clusters
        }


@pytest.fixture()
def synthetic_dataset(tmp_path: Path):
    data = {
        "dataset_name": "synthetic",
        "dataset_context": {"species": "Homo sapiens"},
        "clusters": [
            {"cluster_id": "0", "ground_truth": "B cell", "markers": ["MS4A1"]},
            {"cluster_id": "1", "ground_truth": "T cell", "markers": ["CD3D"]},
        ],
    }
    dataset_path = tmp_path / "dataset.json"
    dataset_path.write_text(json.dumps(data), encoding="utf-8")
    return data, dataset_path


def test_run_gpt_benchmark_metrics(synthetic_dataset):
    get_settings().validation_min_marker_overlap = 1
    data, _ = synthetic_dataset
    annotator = DummyAnnotator()
    result = run_gpt_benchmark(data, annotator, MARKER_DB)
    assert isinstance(result, BenchmarkResult)
    assert result.accuracy == 1.0
    assert result.top3_recall >= result.accuracy
    assert len(result.per_class) > 0


def test_load_and_run(synthetic_dataset):
    get_settings().validation_min_marker_overlap = 1
    _, dataset_path = synthetic_dataset
    annotator = DummyAnnotator()
    result = load_and_run(dataset_path, annotator, MARKER_DB)
    assert result.macro_f1 >= 0.0
    assert result.metadata["runtime_seconds"] >= 0.0


def test_marker_overlap_baseline(synthetic_dataset):
    get_settings().validation_min_marker_overlap = 1
    data, _ = synthetic_dataset
    baseline = run_marker_overlap_baseline(data, MARKER_DB)
    assert baseline.model_name == "marker_overlap_baseline"
    assert baseline.unknown_rate <= 1.0


def test_render_markdown_report_highlights_regression():
    result = {
        "dataset": "synthetic",
        "models": [
            {
                "model_name": "gpt",
                "accuracy": 0.70,
                "macro_f1": 0.65,
                "top3_recall": 0.80,
                "unknown_rate": 0.10,
                "flag_precision": 0.75,
                "time_per_cluster": 0.1,
                "per_class": {},
                "predictions": [],
            }
        ],
    }
    markdown = render_dataset_report(result, deltas={"gpt": {"accuracy": -0.06}})
    assert "**-6.00%**" in markdown


def test_render_sparkline_csv_formats_history():
    csv_text = render_sparkline_csv(
        "synthetic",
        [
            {"date": "20240101", "accuracy": 0.8, "macro_f1": 0.75},
            {"date": "20240201", "accuracy": 0.85, "macro_f1": 0.78},
        ],
    )
    lines = csv_text.strip().splitlines()
    assert lines[0] == "dataset,date,accuracy,macro_f1"
    assert lines[1].startswith("synthetic,20240101,0.8000")
