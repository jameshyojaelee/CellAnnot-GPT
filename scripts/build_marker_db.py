"""CLI entry point for building the GPT Cell Annotator marker knowledge base."""

from __future__ import annotations

import argparse
from pathlib import Path

from backend.data_ingest.marker_loader import (
    ChecksumMismatchError,
    MarkerDataLoader,
    SourceResolutionError,
    default_sources,
    load_sources_from_yaml,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the marker gene knowledge base")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory where parquet/sqlite outputs will be written.",
    )
    parser.add_argument(
        "--skip-parquet",
        action="store_true",
        help="Do not write the parquet artifact.",
    )
    parser.add_argument(
        "--skip-sqlite",
        action="store_true",
        help="Do not write the SQLite artifact.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/marker_sources.yaml"),
        help="YAML file listing marker data sources.",
    )
    parser.add_argument(
        "--local-only",
        action="store_true",
        help="Use only local files; skip network downloads.",
    )
    parser.add_argument(
        "--verify-checksums",
        action="store_true",
        help="Fail if downloaded sources do not match configured checksums.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.config:
        sources = load_sources_from_yaml(args.config)
    else:
        sources = default_sources()

    loader = MarkerDataLoader(
        sources,
        storage_dir=args.output_dir,
        parquet_path=args.output_dir / "marker_db.parquet",
        sqlite_path=args.output_dir / "marker_db.sqlite",
    )

    write_parquet = not args.skip_parquet
    write_sqlite = not args.skip_sqlite

    try:
        df = loader.run(
            write_parquet=write_parquet,
            write_sqlite=write_sqlite,
            local_only=args.local_only,
            enforce_checksums=args.verify_checksums,
        )
    except SourceResolutionError as exc:
        raise SystemExit(
            f"[gpt-cell-annotator] Unable to resolve any marker sources: {exc}"
        ) from exc
    except ChecksumMismatchError as exc:
        raise SystemExit(f"[gpt-cell-annotator] Checksum verification failed: {exc}") from exc

    message = (
        f"Ingested {len(df)} marker records from {len(loader.sources)} sources "
        f"into {args.output_dir}"
    )
    print(message)


if __name__ == "__main__":
    main()
