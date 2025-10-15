"""Console entrypoints for GPT Cell Annotator."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import pandas as pd
from rich.console import Console
from rich.table import Table

from backend.data_ingest.marker_loader import (
    MarkerDataLoader,
    load_sources_from_yaml,
)
from backend.validation.crosscheck import crosscheck_batch
from backend.validation.report import build_structured_report
from backend.llm.annotator import Annotator
from config.settings import get_settings
from gpt_cell_annotator import assets
from gpt_cell_annotator.scanpy import report_to_dataframe

console = Console()


def _package_version() -> str:
    try:
        return version("gpt-cell-annotator")
    except PackageNotFoundError:
        return "0.0.0"


def _prepare_environment(args: argparse.Namespace) -> tuple[Path, Path]:
    """Materialise offline assets and normalise environment variables."""

    home_dir = assets.get_asset_home(args.assets_home)
    assets.ensure_all_assets(home=home_dir)

    data_dir = Path(args.data_dir) if args.data_dir else home_dir / "data/processed"
    assets.ensure_marker_database(home=home_dir, target_dir=data_dir)
    os.environ.setdefault("GPT_CELL_ANNOTATOR_DATA_DIR", str(data_dir))
    if args.offline:
        os.environ.setdefault("OPENAI_API_KEY", "")
    return home_dir, data_dir


def _parse_markers(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return []
        try:
            parsed = json.loads(text)
            if isinstance(parsed, list):
                candidates = parsed
            else:
                candidates = [parsed]
        except json.JSONDecodeError:
            separators = [",", ";", "|"]
            for sep in separators:
                if sep in text:
                    candidates = [part.strip() for part in text.split(sep)]
                    break
            else:
                candidates = [text]
    elif isinstance(value, (list, tuple, set)):
        candidates = list(value)
    else:
        candidates = [value]

    markers: list[str] = []
    for candidate in candidates:
        symbol = str(candidate).strip()
        if symbol:
            markers.append(symbol.upper())
    return markers


def _load_marker_payload(
    path: Path,
    *,
    cluster_column: str,
    markers_column: str,
) -> list[dict[str, Any]]:
    frame = pd.read_csv(path)
    if cluster_column not in frame.columns:
        raise KeyError(f"Missing '{cluster_column}' column in {path}")
    if markers_column not in frame.columns:
        raise KeyError(f"Missing '{markers_column}' column in {path}")

    payload: list[dict[str, Any]] = []
    for row in frame.itertuples():
        cluster_id = getattr(row, cluster_column)
        markers_raw = getattr(row, markers_column)
        markers = _parse_markers(markers_raw)
        payload.append(
            {
                "cluster_id": str(cluster_id),
                "markers": markers,
            }
        )
    return payload


def _print_report(report: dict[str, Any]) -> None:
    table = Table(title="GPT Cell Annotator", expand=True)
    table.add_column("Cluster")
    table.add_column("Primary Label")
    table.add_column("Confidence", justify="center")
    table.add_column("Status", justify="center")
    table.add_column("Warnings")

    for cluster in report.get("clusters", []):
        annotation = cluster.get("annotation") or {}
        warnings = "; ".join(cluster.get("warnings") or [])
        table.add_row(
            str(cluster.get("cluster_id")),
            annotation.get("primary_label") or "Unknown",
            cluster.get("confidence") or "",
            cluster.get("status") or "",
            warnings,
        )

    console.print(table)


def cmd_annotate(args: argparse.Namespace) -> int:
    home_dir, data_dir = _prepare_environment(args)
    assets.ensure_demo_files(home=home_dir)

    try:
        input_path = assets.resolve_path(args.input, home=home_dir)
    except FileNotFoundError as exc:
        console.print(f"[red]Unable to locate marker CSV: {args.input}[/red]")
        raise SystemExit(2) from exc

    payload = _load_marker_payload(
        input_path,
        cluster_column=args.cluster_column,
        markers_column=args.markers_column,
    )

    dataset_context: dict[str, str] = {"species": args.species}
    if args.tissue:
        dataset_context["tissue"] = args.tissue

    settings = get_settings()
    if args.offline:
        annotator = Annotator(settings=settings, force_mock=True)
    else:
        annotator = Annotator(settings=settings)

    raw_result = annotator.annotate_batch(payload, dataset_context)
    marker_db_path = data_dir / "marker_db.parquet"
    marker_db = pd.read_parquet(marker_db_path)

    annotations: list[dict[str, Any]] = []
    for cluster in payload:
        cid = cluster["cluster_id"]
        cluster_result = dict(raw_result.get(cid) or {})
        cluster_result.setdefault("primary_label", "Unknown or Novel")
        cluster_result.setdefault("confidence", "Unknown")
        cluster_result.setdefault("rationale", "")
        annotations.append(
            {
                **cluster_result,
                "cluster_id": cid,
                "markers": cluster["markers"],
            }
        )

    crosschecked = crosscheck_batch(
        annotations,
        marker_db,
        species=args.species,
        tissue=args.tissue,
        min_support=settings.validation_min_marker_overlap,
    )
    report = build_structured_report(annotations, crosschecked).model_dump()
    _print_report(report)

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if args.out_csv:
        df = report_to_dataframe(report)
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out_csv, index=False)

    return 0


def _rewrite_local_paths(sources: list, home_dir: Path) -> None:
    assets.ensure_demo_files(home=home_dir)
    for cfg in sources:
        if getattr(cfg, "local_path", None) and not cfg.local_path.exists():
            try:
                resolved = assets.resolve_path(cfg.local_path.as_posix(), home=home_dir)
            except FileNotFoundError:
                continue
            cfg.local_path = resolved


def cmd_build_db(args: argparse.Namespace) -> int:
    home_dir, default_data_dir = _prepare_environment(args)
    output_dir = Path(args.output_dir) if args.output_dir else default_data_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.config:
        config_path = Path(args.config)
    else:
        config_path = assets.resolve_path("config/marker_sources.yaml", home=home_dir)

    sources = load_sources_from_yaml(config_path)
    _rewrite_local_paths(sources, home_dir)

    loader = MarkerDataLoader(
        sources,
        output_dir,
        parquet_path=output_dir / "marker_db.parquet",
        sqlite_path=output_dir / "marker_db.sqlite",
    )
    try:
        df = loader.run(
            write_parquet=not args.skip_parquet,
            write_sqlite=not args.skip_sqlite,
            local_only=args.offline or args.local_only,
            enforce_checksums=args.verify_checksums,
        )
    finally:
        loader.close()

    console.print(
        f"[green]Ingested {len(df)} marker rows into {output_dir}[/green]",
    )
    return 0


def cmd_scanpy(args: argparse.Namespace) -> int:
    home_dir, data_dir = _prepare_environment(args)
    cmd_args = args.scanpy_args or ["--help"]
    os.environ.setdefault("GPT_CELL_ANNOTATOR_ASSETS_HOME", str(home_dir))
    sys.argv = ["gca scanpy"] + cmd_args
    from gpt_cell_annotator import scanpy

    return scanpy.main(cmd_args)


def cmd_api(args: argparse.Namespace) -> int:
    _prepare_environment(args)
    command = [
        "uvicorn",
        args.app,
        "--host",
        args.host,
        "--port",
        str(args.port),
    ]
    if args.reload:
        command.append("--reload")
    return subprocess.call(command)


def _build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--data-dir",
        type=Path,
        help="Directory for processed marker database artifacts.",
    )
    common.add_argument(
        "--assets-home",
        type=Path,
        help="Directory to cache bundled offline assets.",
    )
    common.add_argument(
        "--offline",
        "--mock",
        dest="offline",
        action="store_true",
        default=False,
        help="Force heuristic mock mode and avoid network access.",
    )

    parser = argparse.ArgumentParser(
        prog="gca",
        description="GPT Cell Annotator command line tools.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=_package_version())

    subparsers = parser.add_subparsers(dest="command", required=True)

    annotate_parser = subparsers.add_parser(
        "annotate",
        parents=[common],
        help="Annotate clusters from a marker CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    annotate_parser.add_argument(
        "input",
        help="Marker CSV path or bundled asset (e.g. data/demo/pbmc_markers.csv).",
    )
    annotate_parser.add_argument(
        "--cluster-column",
        default="cluster_id",
        help="Column containing cluster identifiers.",
    )
    annotate_parser.add_argument(
        "--markers-column",
        default="markers",
        help="Column containing marker gene lists.",
    )
    annotate_parser.add_argument(
        "--species",
        default="Homo sapiens",
        help="Species context for annotation.",
    )
    annotate_parser.add_argument(
        "--tissue",
        help="Optional tissue or compartment context.",
    )
    annotate_parser.add_argument(
        "--out-json",
        type=Path,
        help="Write the structured annotation report to JSON.",
    )
    annotate_parser.add_argument(
        "--out-csv",
        type=Path,
        help="Write a summary CSV of the annotations.",
    )
    annotate_parser.set_defaults(func=cmd_annotate)

    build_parser = subparsers.add_parser(
        "build-db",
        parents=[common],
        help="Build the marker knowledge base artifacts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    build_parser.add_argument(
        "--config",
        type=Path,
        help="Override the marker source configuration file.",
    )
    build_parser.add_argument(
        "--output-dir",
        type=Path,
        help="Destination directory for the generated artifacts.",
    )
    build_parser.add_argument(
        "--skip-parquet",
        action="store_true",
        help="Skip writing the Parquet artifact.",
    )
    build_parser.add_argument(
        "--skip-sqlite",
        action="store_true",
        help="Skip writing the SQLite artifact.",
    )
    build_parser.add_argument(
        "--verify-checksums",
        action="store_true",
        help="Strictly enforce configured checksums for downloads.",
    )
    build_parser.add_argument(
        "--local-only",
        action="store_true",
        help="Restrict ingestion to local files even when URLs are available.",
    )
    build_parser.set_defaults(func=cmd_build_db)

    scanpy_parser = subparsers.add_parser(
        "scanpy",
        parents=[common],
        help="Proxy subcommands for gpt_cell_annotator.scanpy.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    scanpy_parser.add_argument(
        "scanpy_args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to `python -m gpt_cell_annotator.scanpy`.",
    )
    scanpy_parser.set_defaults(func=cmd_scanpy)

    api_parser = subparsers.add_parser(
        "api",
        parents=[common],
        help="Run the FastAPI service.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    api_parser.add_argument(
        "--host",
        default="0.0.0.0",
        help="Interface to bind the API server to.",
    )
    api_parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port to expose the API server on.",
    )
    api_parser.add_argument(
        "--app",
        default="backend.api.main:app",
        help="Import path for the ASGI application.",
    )
    api_parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable autoreload (development mode).",
    )
    api_parser.set_defaults(func=cmd_api)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
