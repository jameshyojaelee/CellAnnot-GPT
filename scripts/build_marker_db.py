"""CLI entry point for building the CellAnnot-GPT marker knowledge base."""

from __future__ import annotations

import argparse
from pathlib import Path

from backend.data_ingest.marker_loader import MarkerDataLoader, default_sources


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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    loader = MarkerDataLoader(
        default_sources(),
        storage_dir=args.output_dir,
        parquet_path=args.output_dir / "marker_db.parquet",
        sqlite_path=args.output_dir / "marker_db.sqlite",
    )

    write_parquet = not args.skip_parquet
    write_sqlite = not args.skip_sqlite

    df = loader.run(write_parquet=write_parquet, write_sqlite=write_sqlite)
    print(  # noqa: T201
        f"Ingested {len(df)} marker records from {len(loader.sources)} sources into {args.output_dir}"
    )


if __name__ == "__main__":
    main()
