"""Benchmark utilities for evaluating GPT Cell Annotator annotations."""

from __future__ import annotations

import json
import time
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, precision_recall_fscore_support

from backend.llm.annotator import Annotator
from backend.validation.crosscheck import crosscheck_batch
from backend.validation.report import build_structured_report
from config.settings import get_settings


@dataclass
class PredictionRecord:
    cluster_id: str
    ground_truth: str
    primary_label: str | None
    alternatives: list[str] = field(default_factory=list)
    status: str = "supported"
    confidence: str | None = None
    flag_reasons: list[str] = field(default_factory=list)


@dataclass
class BenchmarkResult:
    model_name: str
    dataset_name: str
    accuracy: float
    macro_f1: float
    top3_recall: float
    unknown_rate: float
    flag_precision: float
    time_per_cluster: float
    per_class: dict[str, dict[str, float]]
    predictions: list[dict[str, Any]]
    metadata: dict[str, Any] = field(default_factory=dict)


def load_dataset(dataset_path: Path) -> dict[str, Any]:
    with dataset_path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _extract_predictions(
    report: dict[str, Any],
    ground_truth_lookup: dict[str, str],
) -> list[PredictionRecord]:
    records: list[PredictionRecord] = []
    for cluster in report.get("clusters", []):
        annotation = cluster.get("annotation", {})
        alternatives = [
            alt.get("label") for alt in annotation.get("alternatives", []) if isinstance(alt, dict)
        ]
        records.append(
            PredictionRecord(
                cluster_id=str(cluster.get("cluster_id")),
                ground_truth=ground_truth_lookup.get(str(cluster.get("cluster_id")), ""),
                primary_label=annotation.get("primary_label"),
                alternatives=[label for label in alternatives if label],
                status=cluster.get("status", "flagged"),
                confidence=annotation.get("confidence"),
                flag_reasons=cluster.get("warnings", []),
            )
        )
    return records


def _compute_metrics(records: Iterable[PredictionRecord], runtime_seconds: float) -> dict[str, Any]:
    predictions = list(records)
    ground_truth = [rec.ground_truth for rec in predictions]
    predicted = [rec.primary_label or "Unknown or Novel" for rec in predictions]

    accuracy = accuracy_score(ground_truth, predicted)
    macro_f1 = f1_score(ground_truth, predicted, average="macro", zero_division=0)
    precision, recall, f1, labels = precision_recall_fscore_support(
        ground_truth, predicted, zero_division=0
    )
    per_class = {
        str(label): {"precision": float(p), "recall": float(r), "f1": float(f)}
        for label, p, r, f in zip(labels, precision, recall, f1, strict=True)
    }

    top3_hits = 0
    unknown_hits = 0
    flagged_total = 0
    flagged_correct = 0

    for rec in predictions:
        candidates = [label for label in [rec.primary_label, *rec.alternatives] if label]
        candidates = candidates[:3]
        if rec.ground_truth in candidates:
            top3_hits += 1
        if (rec.primary_label or "").lower().startswith("unknown"):
            unknown_hits += 1
        if rec.status != "supported":
            flagged_total += 1
            if rec.primary_label != rec.ground_truth:
                flagged_correct += 1

    top3_recall = top3_hits / len(predictions) if predictions else 0.0
    unknown_rate = unknown_hits / len(predictions) if predictions else 0.0
    flag_precision = flagged_correct / flagged_total if flagged_total else 0.0
    time_per_cluster = runtime_seconds / len(predictions) if predictions else 0.0

    return {
        "accuracy": float(accuracy),
        "macro_f1": float(macro_f1),
        "top3_recall": float(top3_recall),
        "unknown_rate": float(unknown_rate),
        "flag_precision": float(flag_precision),
        "time_per_cluster": float(time_per_cluster),
        "per_class": per_class,
    }


def run_gpt_benchmark(
    dataset: dict[str, Any],
    annotator: Annotator,
    marker_db: pd.DataFrame,
    *,
    model_name: str = "gpt_cell_annotator",
) -> BenchmarkResult:
    clusters = dataset["clusters"]
    context = dataset.get("dataset_context")
    payload = [
        {"cluster_id": str(cluster["cluster_id"]), "markers": cluster["markers"]}
        for cluster in clusters
    ]

    start = time.perf_counter()
    batch_result = annotator.annotate_batch(payload, context)
    runtime = time.perf_counter() - start

    annotations = []
    ground_truth_lookup = {
        str(cluster["cluster_id"]): cluster["ground_truth"] for cluster in clusters
    }
    for cluster in clusters:
        cluster_id = str(cluster["cluster_id"])
        annotations.append(
            {
                **(batch_result.get(cluster_id) or {}),
                "cluster_id": cluster_id,
                "markers": cluster["markers"],
            }
        )

    settings = get_settings()
    crosschecked = crosscheck_batch(
        annotations,
        marker_db,
        species=(context or {}).get("species"),
        tissue=(context or {}).get("tissue"),
        min_support=settings.validation_min_marker_overlap,
    )
    report = build_structured_report(annotations, crosschecked).model_dump()
    records = _extract_predictions(report, ground_truth_lookup)
    metrics = _compute_metrics(records, runtime)

    return BenchmarkResult(
        model_name=model_name,
        dataset_name=dataset.get("dataset_name", "unknown"),
        accuracy=metrics["accuracy"],
        macro_f1=metrics["macro_f1"],
        top3_recall=metrics["top3_recall"],
        unknown_rate=metrics["unknown_rate"],
        flag_precision=metrics["flag_precision"],
        time_per_cluster=metrics["time_per_cluster"],
        per_class=metrics["per_class"],
        predictions=[record.__dict__ for record in records],
        metadata={"runtime_seconds": runtime},
    )


def run_marker_overlap_baseline(
    dataset: dict[str, Any],
    marker_db: pd.DataFrame,
    *,
    model_name: str = "marker_overlap_baseline",
) -> BenchmarkResult:
    context = dataset.get("dataset_context") or {}
    markers_grouped = marker_db.groupby("cell_type")
    cluster_predictions: list[PredictionRecord] = []
    start = time.perf_counter()

    for cluster in dataset["clusters"]:
        markers = {marker.upper() for marker in cluster["markers"]}
        best_label = "Unknown or Novel"
        best_overlap = 0
        overlaps: list[tuple[str, int]] = []
        for label, df in markers_grouped:
            df_species = df
            if context.get("species"):
                df_species = df_species[
                    df_species["species"].str.lower() == context["species"].lower()
                ]
            if df_species.empty:
                continue
            overlap = len(markers & set(df_species["gene_symbol"].str.upper()))
            if overlap:
                overlaps.append((label, overlap))
                if overlap > best_overlap:
                    best_label = label
                    best_overlap = overlap

        overlaps.sort(key=lambda item: item[1], reverse=True)
        alternatives = [label for label, _ in overlaps[1:3]]
        if best_overlap >= 3:
            confidence = "High"
        elif best_overlap == 2:
            confidence = "Medium"
        else:
            confidence = "Low"

        cluster_predictions.append(
            PredictionRecord(
                cluster_id=str(cluster["cluster_id"]),
                ground_truth=cluster["ground_truth"],
                primary_label=best_label if best_overlap else "Unknown or Novel",
                alternatives=alternatives,
                status="supported" if best_overlap else "flagged",
                confidence=confidence,
                flag_reasons=[] if best_overlap else ["No marker overlap with knowledge base"],
            )
        )

    runtime = time.perf_counter() - start
    metrics = _compute_metrics(cluster_predictions, runtime)
    return BenchmarkResult(
        model_name=model_name,
        dataset_name=dataset.get("dataset_name", "unknown"),
        accuracy=metrics["accuracy"],
        macro_f1=metrics["macro_f1"],
        top3_recall=metrics["top3_recall"],
        unknown_rate=metrics["unknown_rate"],
        flag_precision=metrics["flag_precision"],
        time_per_cluster=metrics["time_per_cluster"],
        per_class=metrics["per_class"],
        predictions=[record.__dict__ for record in cluster_predictions],
        metadata={"runtime_seconds": runtime},
    )


def load_and_run(
    dataset_path: Path,
    annotator: Annotator,
    marker_db: pd.DataFrame,
) -> BenchmarkResult:
    dataset = load_dataset(dataset_path)
    return run_gpt_benchmark(dataset, annotator, marker_db)


__all__ = [
    "BenchmarkResult",
    "PredictionRecord",
    "load_and_run",
    "load_dataset",
    "run_gpt_benchmark",
    "run_marker_overlap_baseline",
]
