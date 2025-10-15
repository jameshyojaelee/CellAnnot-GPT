"""Plotting utilities for benchmark reporting."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

CONFIDENCE_MAP = {"High": 0.9, "Medium": 0.6, "Low": 0.3, "Unknown": 0.0}


def _ensure_output(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def plot_confusion_matrix(predictions: Iterable[dict], output_path: Path) -> None:
    truth_labels = sorted({row["ground_truth"] for row in predictions})
    pred_labels = sorted({row["primary_label"] for row in predictions})
    labels = sorted(set(truth_labels) | set(pred_labels))
    matrix = np.zeros((len(labels), len(labels)), dtype=int)

    label_to_idx = {label: idx for idx, label in enumerate(labels)}
    for row in predictions:
        truth_idx = label_to_idx[row["ground_truth"]]
        pred_idx = label_to_idx[row["primary_label"]]
        matrix[truth_idx, pred_idx] += 1

    plt.figure(figsize=(6, 5))
    sns.heatmap(matrix, annot=True, fmt="d", cmap="Blues", xticklabels=labels, yticklabels=labels)
    plt.xlabel("Predicted")
    plt.ylabel("Ground truth")
    plt.title("Confusion matrix")
    _ensure_output(output_path)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_precision_recall_bars(per_class: dict[str, dict[str, float]], output_path: Path) -> None:
    labels = list(per_class.keys())
    precision = [per_class[label]["precision"] for label in labels]
    recall = [per_class[label]["recall"] for label in labels]

    x = np.arange(len(labels))
    width = 0.35

    plt.figure(figsize=(max(6, len(labels) * 1.2), 4))
    plt.bar(x - width / 2, precision, width, label="Precision")
    plt.bar(x + width / 2, recall, width, label="Recall")
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.ylim(0, 1)
    plt.ylabel("Score")
    plt.title("Per-class precision & recall")
    plt.legend()
    _ensure_output(output_path)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_calibration(predictions: Iterable[dict], output_path: Path) -> None:
    scores = []
    outcomes = []
    for row in predictions:
        confidence = row.get("confidence") or "Unknown"
        score = CONFIDENCE_MAP.get(confidence.title(), 0.0)
        scores.append(score)
        outcomes.append(1.0 if row.get("primary_label") == row.get("ground_truth") else 0.0)

    if not scores:
        return

    plt.figure(figsize=(5, 4))
    plt.scatter(scores, outcomes, alpha=0.6)
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
    plt.xlabel("Confidence (mapped)")
    plt.ylabel("Empirical accuracy")
    plt.title("Calibration curve (coarse)")
    _ensure_output(output_path)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


__all__ = [
    "plot_calibration",
    "plot_confusion_matrix",
    "plot_precision_recall_bars",
]
