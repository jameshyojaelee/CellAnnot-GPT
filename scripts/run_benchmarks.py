#!/usr/bin/env python3
"""Run GPT Cell Annotator benchmarks against curated marker datasets."""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd

from backend.llm.annotator import Annotator
from backend.llm.annotator import MockAnnotator
from config.settings import Settings
from evaluation.benchmark_runner import (
    BenchmarkResult,
    load_dataset,
    run_gpt_benchmark,
    run_marker_overlap_baseline,
)
from evaluation.plots import (
    plot_calibration,
    plot_confusion_matrix,
    plot_precision_recall_bars,
)
from evaluation.report_templates import render_dataset_report

logger = logging.getLogger("gpt_cell_annotator.benchmarks")


def _discover_datasets(dataset_dir: Path, filters: set[str] | None) -> list[Path]:
    datasets = sorted(dataset_dir.glob("*.json"))
    if filters:
        datasets = [
            path for path in datasets if path.stem in filters or path.name in filters
        ]
    return datasets


def _load_marker_db(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Marker DB not found at {path}. Run scripts/build_marker_db.py first."
        )
    df = pd.read_parquet(path)
    required_columns = {
        "cell_type",
        "gene_symbol",
        "species",
        "tissue",
        "ontology_id",
    }
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Marker DB missing columns: {', '.join(sorted(missing))}")
    return df


def _result_to_dict(result: BenchmarkResult) -> dict:
    data = asdict(result)
    data["per_class"] = {label: metrics for label, metrics in result.per_class.items()}
    return data


def _write_model_outputs(
    dataset_dir: Path,
    result: BenchmarkResult,
) -> None:
    model_dir = dataset_dir / result.model_name
    model_dir.mkdir(parents=True, exist_ok=True)
    json_path = model_dir / "metrics.json"
    json_path.write_text(json.dumps(_result_to_dict(result), indent=2), encoding="utf-8")

    predictions_path = model_dir / "predictions.json"
    predictions_path.write_text(json.dumps(result.predictions, indent=2), encoding="utf-8")

    predictions = result.predictions
    if predictions:
        plot_confusion_matrix(predictions, model_dir / "confusion_matrix.png")
        plot_precision_recall_bars(result.per_class, model_dir / "precision_recall.png")
        plot_calibration(predictions, model_dir / "calibration.png")


def _build_dataset_summary(dataset_name: str, model_results: Iterable[BenchmarkResult]) -> dict:
    return {
        "dataset": dataset_name,
        "models": [
            {
                **_result_to_dict(result),
            }
            for result in model_results
        ],
    }


def run_benchmarks(
    datasets: Iterable[Path],
    *,
    output_dir: Path,
    marker_db: pd.DataFrame,
    use_mock: bool,
    baselines: list[str],
) -> list[dict]:
    settings = Settings()
    if use_mock:
        settings.openai_api_key = ""
    annotator = Annotator(settings=settings, mock_backend=MockAnnotator())

    run_id = datetime.utcnow().strftime("%Y%m%d-%H%M%S")
    out_root = output_dir / run_id
    out_root.mkdir(parents=True, exist_ok=True)

    dataset_summaries: list[dict] = []

    for dataset_path in datasets:
        dataset = load_dataset(dataset_path)
        dataset_name = dataset.get("dataset_name", dataset_path.stem)
        logger.info("Running dataset %s", dataset_name)

        dataset_dir = out_root / dataset_name
        dataset_dir.mkdir(parents=True, exist_ok=True)

        model_results: list[BenchmarkResult] = []
        gpt_result = run_gpt_benchmark(dataset, annotator, marker_db)
        model_results.append(gpt_result)
        _write_model_outputs(dataset_dir, gpt_result)

        if "marker_overlap" in baselines:
            baseline = run_marker_overlap_baseline(dataset, marker_db)
            model_results.append(baseline)
            _write_model_outputs(dataset_dir, baseline)

        dataset_summary = _build_dataset_summary(dataset_name, model_results)
        dataset_summaries.append(dataset_summary)

        markdown = render_dataset_report(dataset_summary)
        (dataset_dir / "README.md").write_text(markdown, encoding="utf-8")

    summary_path = out_root / "summary.json"
    summary_path.write_text(json.dumps(dataset_summaries, indent=2), encoding="utf-8")

    history_csv = out_root / "history.csv"
    history_rows = []
    for summary in dataset_summaries:
        for model in summary["models"]:
            history_rows.append(
                {
                    "dataset": summary["dataset"],
                    "model": model["model_name"],
                    "date": run_id,
                    "accuracy": model["accuracy"],
                    "macro_f1": model["macro_f1"],
                }
            )
    pd.DataFrame(history_rows).to_csv(history_csv, index=False)

    logger.info("Reports written to %s", out_root)
    return dataset_summaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark GPT Cell Annotator.")
    parser.add_argument(
        "--datasets",
        help="Comma-separated dataset names (defaults to all JSON files in evaluation/datasets/).",
    )
    parser.add_argument(
        "--dataset-dir",
        default="evaluation/datasets",
        type=Path,
        help="Directory containing dataset metadata JSON files.",
    )
    parser.add_argument(
        "--output",
        default="docs/reports",
        type=Path,
        help="Directory where benchmark reports will be written.",
    )
    parser.add_argument(
        "--marker-db",
        default="data/processed/marker_db.parquet",
        type=Path,
        help="Path to marker knowledge base parquet artifact.",
    )
    parser.add_argument(
        "--baselines",
        default="marker_overlap",
        help="Comma-separated baselines to run (marker_overlap).",
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Force mock annotator (offline mode).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING, ERROR).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))

    filters = set(filter(None, (args.datasets or "").split(","))) if args.datasets else None
    dataset_dir: Path = args.dataset_dir
    dataset_paths = _discover_datasets(dataset_dir, filters)
    if not dataset_paths:
        logger.error("No datasets found in %s", dataset_dir)
        return 1

    try:
        marker_db = _load_marker_db(args.marker_db)
    except (FileNotFoundError, ValueError) as exc:
        logger.error("Marker DB error: %s", exc)
        return 1

    baselines = [item.strip() for item in args.baselines.split(",") if item.strip()]
    run_benchmarks(
        dataset_paths,
        output_dir=args.output,
        marker_db=marker_db,
        use_mock=args.mock,
        baselines=baselines,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
