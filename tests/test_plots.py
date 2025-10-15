from __future__ import annotations

from pathlib import Path

from evaluation.plots import (
    plot_calibration,
    plot_confusion_matrix,
    plot_precision_recall_bars,
)


def test_plot_helpers(tmp_path: Path) -> None:
    predictions = [
        {
            "cluster_id": "0",
            "ground_truth": "B cell",
            "primary_label": "B cell",
            "confidence": "High",
        },
        {
            "cluster_id": "1",
            "ground_truth": "T cell",
            "primary_label": "Monocyte",
            "confidence": "Low",
        },
    ]
    per_class = {
        "B cell": {"precision": 1.0, "recall": 1.0, "f1": 1.0},
        "T cell": {"precision": 0.0, "recall": 0.0, "f1": 0.0},
    }

    plot_confusion_matrix(predictions, tmp_path / "confusion.png")
    plot_precision_recall_bars(per_class, tmp_path / "pr.png")
    plot_calibration(predictions, tmp_path / "calibration.png")

    assert (tmp_path / "confusion.png").exists()
    assert (tmp_path / "pr.png").exists()
    assert (tmp_path / "calibration.png").exists()
