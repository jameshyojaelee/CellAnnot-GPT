"""Run benchmark datasets, track deltas, and persist reports."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from backend.llm.annotator import Annotator
from config.settings import get_settings
from evaluation.benchmark_runner import load_and_run
from evaluation.report_templates import (
    render_markdown_report,
    render_sparkline_csv,
    render_text_confusion_matrix,
)

BENCHMARK_REGRESSION_THRESHOLD = 0.05  # 5 percentage points


def discover_datasets(directory: Path, patterns: Optional[Iterable[str]] = None) -> List[Path]:
    if not patterns:
        return sorted(p for p in directory.glob("*.json") if p.is_file())

    resolved: Dict[Path, None] = {}
    for pattern in patterns:
        candidate = Path(pattern)
        matches: List[Path] = []
        if candidate.is_absolute():
            matches = sorted(p for p in candidate.parent.glob(candidate.name) if p.is_file())
        else:
            direct = directory / candidate
            if direct.exists():
                matches = [direct]
            else:
                matches = sorted(p for p in directory.glob(pattern) if p.is_file())
        for match in matches:
            resolved[match.resolve()] = None
    return sorted(Path(path) for path in resolved.keys())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CellAnnot-GPT benchmarks")
    parser.add_argument(
        "--datasets-dir",
        type=Path,
        default=Path("evaluation/datasets"),
        help="Directory containing benchmark JSON configs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/reports"),
        help="Directory where reports will be stored.",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        help="Specific dataset filenames or glob patterns (defaults to all *.json in directory).",
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Force annotator into mock mode by clearing OPENAI_API_KEY.",
    )
    return parser.parse_args()


def load_previous_metrics(previous_dir: Path, dataset_name: str) -> Optional[Dict[str, Optional[float]]]:
    previous_file = previous_dir / f"{dataset_name}.json"
    if not previous_file.exists():
        return None
    try:
        with previous_file.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
        previous_accuracy = data.get("accuracy")
        previous_macro = data.get("macro_f1")
        return {
            "accuracy": float(previous_accuracy) if previous_accuracy is not None else None,
            "macro_f1": float(previous_macro) if previous_macro is not None else None,
            "run_date": data.get("run_date"),
        }
    except json.JSONDecodeError:
        return None


def write_outputs(
    output_prefix: Path,
    dataset_name: str,
    payload: Dict[str, Any],
    deltas: Dict[str, Optional[float]],
    history: List[Dict[str, float]],
) -> None:
    output_prefix.with_suffix(".json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    markdown_payload = render_markdown_report(payload, deltas=deltas, regression_threshold=BENCHMARK_REGRESSION_THRESHOLD)
    markdown_payload += "\n\n" + render_text_confusion_matrix(payload)
    output_prefix.with_suffix(".md").write_text(markdown_payload, encoding="utf-8")

    sparkline_csv = render_sparkline_csv(dataset_name, history)
    output_prefix.with_suffix(".csv").write_text(sparkline_csv, encoding="utf-8")


def update_latest_dir(latest_dir: Path, timestamp_dir: Path, dataset_names: Iterable[str], summary_files: Iterable[Path]) -> None:
    if latest_dir.exists():
        shutil.rmtree(latest_dir)
    latest_dir.mkdir(parents=True, exist_ok=True)

    for dataset in dataset_names:
        for suffix in (".json", ".md", ".csv"):
            src = (timestamp_dir / dataset).with_suffix(suffix)
            if src.exists():
                shutil.copy2(src, latest_dir / f"{dataset}{suffix}")

    for summary in summary_files:
        if summary.exists():
            shutil.copy2(summary, latest_dir / summary.name)


def run() -> None:
    args = parse_args()
    datasets = discover_datasets(args.datasets_dir, args.datasets)
    if not datasets:
        raise SystemExit("No benchmark datasets found.")

    settings = get_settings()
    if args.mock:
        settings.openai_api_key = ""
    annotator = Annotator(settings=settings)

    timestamp = datetime.utcnow().strftime("%Y%m%d")
    base_output = args.output_dir / timestamp
    base_output.mkdir(parents=True, exist_ok=True)

    latest_dir = args.output_dir / "latest"
    diff_summary: Dict[str, Dict[str, Optional[float]]] = {}
    regressions: List[str] = []

    for dataset_path in datasets:
        dataset_name = dataset_path.stem
        result = load_and_run(dataset_path, annotator)

        previous = load_previous_metrics(latest_dir, dataset_name)
        prev_accuracy = previous.get("accuracy") if previous else None
        prev_macro = previous.get("macro_f1") if previous else None

        if prev_accuracy is not None:
            accuracy_delta = result.accuracy - prev_accuracy
        else:
            accuracy_delta = None
        if prev_macro is not None:
            macro_f1_delta = result.macro_f1 - prev_macro
        else:
            macro_f1_delta = None

        json_payload = {
            "dataset": dataset_name,
            "run_date": timestamp,
            "accuracy": result.accuracy,
            "macro_f1": result.macro_f1,
            "per_class": result.per_class,
            "predictions": result.predictions,
            "baseline": previous,
            "delta": {
                "accuracy": accuracy_delta,
                "macro_f1": macro_f1_delta,
            },
        }
        output_prefix = base_output / dataset_name

        history = []
        if prev_accuracy is not None and prev_macro is not None:
            history.append(
                {
                    "date": previous.get("run_date", "previous") if previous else "previous",
                    "accuracy": prev_accuracy,
                    "macro_f1": prev_macro,
                }
            )
        history.append(
            {
                "date": timestamp,
                "accuracy": result.accuracy,
                "macro_f1": result.macro_f1,
            }
        )

        write_outputs(
            output_prefix,
            dataset_name,
            json_payload,
            {"accuracy": accuracy_delta, "macro_f1": macro_f1_delta},
            history,
        )

        diff_summary[dataset_name] = {
            "accuracy": result.accuracy,
            "macro_f1": result.macro_f1,
            "accuracy_delta": accuracy_delta,
            "macro_f1_delta": macro_f1_delta,
        }

        if accuracy_delta is not None and accuracy_delta <= -BENCHMARK_REGRESSION_THRESHOLD:
            regressions.append(dataset_name)

        accuracy_msg = f"{result.accuracy:.2%}"
        if accuracy_delta is not None:
            accuracy_msg += f" (Δ {accuracy_delta:+.2%})"

        macro_msg = f"{result.macro_f1:.2%}"
        if macro_f1_delta is not None:
            macro_msg += f" (Δ {macro_f1_delta:+.2%})"

        print(f"[benchmarks] {dataset_name}: accuracy {accuracy_msg}, macro_f1 {macro_msg}")

    diff_summary_path = base_output / "diff_summary.json"
    diff_summary_path.write_text(json.dumps(diff_summary, indent=2), encoding="utf-8")

    summary_paths = [diff_summary_path]
    dataset_names = [path.stem for path in datasets]
    update_latest_dir(latest_dir, base_output, dataset_names, summary_paths)

    allow_regressions = os.getenv("ALLOW_BENCHMARK_REGRESSIONS", "").lower() == "true"
    if regressions and not allow_regressions:
        raise SystemExit(f"Accuracy regressions detected (>5% drop): {', '.join(regressions)}")


if __name__ == "__main__":
    run()
