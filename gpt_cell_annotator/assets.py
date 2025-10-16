"""Utilities for working with the packaged demo/offline assets."""

from __future__ import annotations

import os
import shutil
from importlib import resources
from importlib.abc import Traversable
from pathlib import Path
from typing import Iterable, Iterator

ASSET_PACKAGE = "gpt_cell_annotator._assets"
HOME_ENV_VAR = "GPT_CELL_ANNOTATOR_HOME"
DEFAULT_HOME_DIRNAME = ".cache/gpt-cell-annotator"

BUNDLED_FILES: tuple[str, ...] = (
    "config/marker_sources.yaml",
    "data/demo/cellmarker_demo.csv",
    "data/demo/curated_demo.json",
    "data/demo/panglaodb_demo.csv",
    "data/demo/pbmc_markers.csv",
    "data/processed/marker_db.parquet",
    "data/processed/marker_db.sqlite",
)

MARKER_DB_FILES: tuple[str, ...] = (
    "data/processed/marker_db.parquet",
    "data/processed/marker_db.sqlite",
)

DEMO_FILES: tuple[str, ...] = (
    "data/demo/panglaodb_demo.csv",
    "data/demo/cellmarker_demo.csv",
    "data/demo/curated_demo.json",
    "data/demo/pbmc_markers.csv",
)


class AssetNotFoundError(FileNotFoundError):
    """Raised when a requested packaged asset is missing."""


def get_asset_home(root: Path | None = None) -> Path:
    """Return the directory where packaged assets should be materialised."""

    if root is not None:
        home = Path(root)
    else:
        override = os.environ.get(HOME_ENV_VAR)
        if override:
            home = Path(override)
        else:
            home = Path.home() / DEFAULT_HOME_DIRNAME
    home.mkdir(parents=True, exist_ok=True)
    return home


def _iter_traversable(path: Traversable) -> Iterator[tuple[str, Traversable]]:
    if path.is_file():
        yield path.name, path
        return
    for child in path.iterdir():
        if child.is_dir():
            for descendant_name, descendant in _iter_traversable(child):
                yield f"{child.name}/{descendant_name}", descendant
        else:
            yield child.name, child


def _copy_traversable(src: Traversable, dest: Path, *, overwrite: bool) -> None:
    if src.is_dir():
        for name, descendant in _iter_traversable(src):
            target = dest / src.name / name if src.name else dest / name
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() and not overwrite:
                continue
            target.write_bytes(descendant.read_bytes())
    else:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and not overwrite:
            return
        dest.write_bytes(src.read_bytes())


def install_asset(relative_path: str | Path, *, home: Path | None = None, overwrite: bool = False) -> Path:
    """Ensure ``relative_path`` is available under the asset home directory."""

    rel = Path(relative_path)
    if rel.is_absolute():
        raise ValueError("relative_path must be relative")

    home_dir = get_asset_home(home)
    destination = home_dir / rel
    if destination.exists() and not overwrite:
        return destination

    try:
        asset = resources.files(ASSET_PACKAGE).joinpath(rel.as_posix())
    except FileNotFoundError as exc:  # pragma: no cover - defensive
        raise AssetNotFoundError(rel) from exc

    if not (asset.is_file() or asset.is_dir()):
        raise AssetNotFoundError(rel)

    if asset.is_dir():
        _copy_traversable(asset, destination.parent, overwrite=overwrite)
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(asset.read_bytes())
    return destination


def ensure_assets(relative_paths: Iterable[str], *, home: Path | None = None) -> Path:
    """Materialise a collection of assets and return the home directory."""

    home_dir = get_asset_home(home)
    for rel in relative_paths:
        install_asset(rel, home=home_dir)
    return home_dir


def ensure_all_assets(*, home: Path | None = None) -> Path:
    """Ensure every bundled asset is materialised."""

    return ensure_assets(BUNDLED_FILES, home=home)


def ensure_marker_database(
    *,
    home: Path | None = None,
    target_dir: Path | None = None,
) -> Path:
    """Ensure the offline marker database artifacts are available."""

    home_dir = ensure_assets(MARKER_DB_FILES, home=home)
    packaged_dir = home_dir / "data/processed"
    if target_dir is None:
        return packaged_dir

    target = Path(target_dir)
    target.mkdir(parents=True, exist_ok=True)
    for rel in MARKER_DB_FILES:
        source = home_dir / rel
        destination = target / Path(rel).name
        if destination.exists():
            try:
                if source.resolve(strict=True) == destination.resolve(strict=True):
                    continue
            except FileNotFoundError:
                pass
        shutil.copy2(source, destination)
    return target


def ensure_demo_files(
    *,
    home: Path | None = None,
    target_dir: Path | None = None,
) -> Path:
    """Ensure demo CSV/JSON assets are available."""

    home_dir = ensure_assets(DEMO_FILES, home=home)
    packaged_dir = home_dir / "data/demo"
    if target_dir is None:
        return packaged_dir

    target = Path(target_dir)
    target.mkdir(parents=True, exist_ok=True)
    for rel in DEMO_FILES:
        source = home_dir / rel
        destination = target / Path(rel).name
        if destination.exists():
            try:
                if source.resolve(strict=True) == destination.resolve(strict=True):
                    continue
            except FileNotFoundError:
                pass
        shutil.copy2(source, destination)
    return target


def resolve_path(path: str | Path, *, home: Path | None = None, ensure: bool = True) -> Path:
    """Resolve ``path`` to a concrete ``Path``, falling back to packaged assets."""

    candidate = Path(path)
    if candidate.expanduser().exists():
        return candidate.expanduser()

    rel = candidate.as_posix().lstrip("./")
    home_dir = get_asset_home(home)
    if ensure:
        try:
            install_asset(rel, home=home_dir)
        except AssetNotFoundError:
            pass
    fallback = home_dir / rel
    if fallback.exists():
        return fallback
    raise FileNotFoundError(f"Unable to resolve asset path: {path}")


__all__ = [
    "AssetNotFoundError",
    "ensure_all_assets",
    "ensure_assets",
    "ensure_demo_files",
    "ensure_marker_database",
    "get_asset_home",
    "install_asset",
    "resolve_path",
]
