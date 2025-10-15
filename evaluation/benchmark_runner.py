"""Benchmark utilities for evaluating CellAnnot-GPT annotations."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.metrics import accuracy_score, f1_score, precision_recall_fscore_support

from backend.llm.annotator import Annotator


@dataclass
class BenchmarkResult:
    accuracy: float
    macro_f1: float
    per_class: dict[str, dict[str, float]]
    predictions: list[dict[str, Any]]


def load_dataset(dataset_path: Path) -> dict[str, Any]:
    with dataset_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def run_benchmark(
    dataset: dict[str, Any],
    annotator: Annotator,
) -> BenchmarkResult:
    clusters = dataset["clusters"]
    predicted_labels = []
    ground_truth_labels = []
    predictions = []

    for cluster in clusters:
        payload = {"cluster_id": cluster["cluster_id"], "markers": cluster["markers"]}
        response = annotator.annotate_cluster(payload, dataset.get("dataset_context"))
        primary_label = response.get("primary_label")
        predicted_labels.append(primary_label)
        ground_truth_labels.append(cluster["ground_truth"])
        predictions.append(
            {
                "cluster_id": cluster["cluster_id"],
                "ground_truth": cluster["ground_truth"],
                "predicted": primary_label,
                "rationale": response.get("rationale"),
            }
        )

    accuracy = accuracy_score(ground_truth_labels, predicted_labels)
    macro_f1 = f1_score(ground_truth_labels, predicted_labels, average="macro", zero_division=0)
    precision, recall, f1, _ = precision_recall_fscore_support(
        ground_truth_labels, predicted_labels, zero_division=0
    )

    per_class_metrics = {}
    classes = np.unique(ground_truth_labels)
    for idx, cls in enumerate(classes):
        per_class_metrics[cls] = {
            "precision": float(precision[idx]),
            "recall": float(recall[idx]),
            "f1": float(f1[idx]),
        }

    return BenchmarkResult(
        accuracy=float(accuracy),
        macro_f1=float(macro_f1),
        per_class=per_class_metrics,
        predictions=predictions,
    )


def load_and_run(dataset_path: Path, annotator: Annotator) -> BenchmarkResult:
    dataset = load_dataset(dataset_path)
    return run_benchmark(dataset, annotator)


__all__ = ["BenchmarkResult", "load_dataset", "run_benchmark", "load_and_run"]
