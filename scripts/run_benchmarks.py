"""Run benchmark datasets and persist reports."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import List

from backend.llm.annotator import Annotator
from evaluation.benchmark_runner import load_and_run
from evaluation.report_templates import render_markdown_report, render_text_confusion_matrix
from config.settings import get_settings


def discover_datasets(directory: Path) -> List[Path]:
    return sorted(p for p in directory.glob("*.json") if p.is_file())


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
        help="Specific dataset filenames to run (defaults to all in directory).",
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Force annotator into mock mode by clearing OPENAI_API_KEY.",
    )
    return parser.parse_args()


def run() -> None:
    args = parse_args()
    datasets = discover_datasets(args.datasets_dir) if not args.datasets else [args.datasets_dir / d for d in args.datasets]

    if not datasets:
        raise SystemExit("No benchmark datasets found.")

    settings = get_settings()
    if args.mock:
        settings.openai_api_key = ""
    annotator = Annotator(settings=settings)

    timestamp = datetime.utcnow().strftime("%Y%m%d")
    base_output = args.output_dir / timestamp
    base_output.mkdir(parents=True, exist_ok=True)

    for dataset_path in datasets:
        result = load_and_run(dataset_path, annotator)
        dataset_name = dataset_path.stem
        output_prefix = base_output / dataset_name

        json_payload = {
            "dataset": dataset_name,
            "accuracy": result.accuracy,
            "macro_f1": result.macro_f1,
            "per_class": result.per_class,
            "predictions": result.predictions,
        }
        output_prefix.with_suffix(".json").write_text(
            json.dumps(json_payload, indent=2),
            encoding="utf-8",
        )

        markdown_payload = render_markdown_report(json_payload)
        markdown_payload += "\n\n" + render_text_confusion_matrix(json_payload)
        output_prefix.with_suffix(".md").write_text(markdown_payload, encoding="utf-8")

        print(f"Saved benchmark report for {dataset_name} -> {output_prefix}")


if __name__ == "__main__":
    run()
