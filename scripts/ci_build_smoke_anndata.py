#!/usr/bin/env python3
"""Generate a minimal AnnData file for CI smoke tests."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def build_dataset(output_dir: Path) -> None:
    """Create a tiny AnnData object and persist it to ``output_dir``."""

    output_dir.mkdir(parents=True, exist_ok=True)
    matrix = np.zeros((6, 4), dtype=float)
    matrix[:3, 0] = 5
    matrix[3:, 1] = 5
    obs = pd.DataFrame({"cluster": ["0", "0", "0", "1", "1", "1"]})
    var = pd.DataFrame(index=["MS4A1", "CD3E", "GNLY", "LYZ"])
    adata = sc.AnnData(matrix, obs=obs, var=var)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, groupby="cluster", n_genes=3)
    adata.write_h5ad(output_dir / "smoke.h5ad")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        msg = "Usage: ci_build_smoke_anndata.py OUTPUT_DIR"
        raise SystemExit(msg)
    build_dataset(Path(argv[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
