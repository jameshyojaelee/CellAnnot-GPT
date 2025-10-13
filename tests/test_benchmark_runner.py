from __future__ import annotations


import pytest

from pathlib import Path
import json

from evaluation.benchmark_runner import BenchmarkResult, load_and_run, run_benchmark


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
