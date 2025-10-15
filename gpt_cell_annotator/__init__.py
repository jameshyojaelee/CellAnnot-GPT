"""Public helpers for integrating GPT Cell Annotator into Python workflows.

Importing heavy submodules (e.g., Scanpy integration) is deferred via lazy
attribute access to avoid circular imports during backend initialisation.
"""

from __future__ import annotations

from types import ModuleType
from typing import Any

__all__ = ["annotate_anndata", "report_to_dataframe"]


def _load_scanpy_module() -> ModuleType:
    from . import scanpy  # local import to avoid eager dependency loading

    return scanpy


def __getattr__(name: str) -> Any:
    if name in __all__:
        module = _load_scanpy_module()
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
