from __future__ import annotations


import pytest

from pathlib import Path
import json

from evaluation.benchmark_runner import BenchmarkResult, load_and_run, run_benchmark
from evaluation.report_templates import render_markdown_report, render_sparkline_csv


class DummyAnnotator:
    def annotate_cluster(self, payload, dataset_context=None):  # type: ignore[override]
        cluster_id = payload["cluster_id"]
        if cluster_id == "0":
            return {"primary_label": "B cell", "rationale": "MS4A1 detected"}
        return {"primary_label": "Monocyte", "rationale": "Default"}


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


def test_run_benchmark_metrics(synthetic_dataset):
    data, _ = synthetic_dataset
    annotator = DummyAnnotator()
    result = run_benchmark(data, annotator)
    assert isinstance(result, BenchmarkResult)
    assert result.accuracy == 0.5
    assert "B cell" in result.per_class


def test_load_and_run(synthetic_dataset):
    _, dataset_path = synthetic_dataset
    annotator = DummyAnnotator()
    result = load_and_run(dataset_path, annotator)
    assert result.macro_f1 >= 0.0


def test_render_markdown_report_highlights_regression():
    result = {
        "accuracy": 0.70,
        "macro_f1": 0.65,
        "per_class": {},
        "predictions": [],
    }
    markdown = render_markdown_report(result, deltas={"accuracy": -0.06, "macro_f1": 0.01})
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
