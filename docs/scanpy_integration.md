# Scanpy Integration Guide

This guide shows how to use GPT Cell Annotator directly from Scanpy notebooks and command-line workflows.

## 1. Prerequisites
- Build the marker knowledge base: `python scripts/build_marker_db.py`
- Install Scanpy support: `poetry install` (or `pip install -r requirements.txt`) — the PyPI package includes `scanpy` as an optional dependency.
- Ensure your dataset has cluster assignments stored in `adata.obs`.

## 2. Annotate within a notebook
```python
import scanpy as sc
from gpt_cell_annotator import annotate_anndata

# Load AnnData object with clustering stored in adata.obs["leiden"]
adata = sc.read_h5ad("data/demo/pbmc_demo.h5ad")

# Optional: compute differential expression if not done already
sc.tl.rank_genes_groups(adata, groupby="leiden", n_genes=5, method="wilcoxon")

# Annotate clusters
adata, report = annotate_anndata(
    adata,
    "leiden",
    species="Homo sapiens",
    tissue="Blood",
)

# Inspect results written to adata.obs
adata.obs[["gptca_label", "gptca_confidence", "gptca_status"]].head()
```

The helper automatically:
1. Computes `rank_genes_groups` if missing (requires Scanpy).
2. Extracts the top N markers per cluster and calls the batch annotator.
3. Cross-checks predictions against the marker DB and writes columns:
   - `gptca_label`
   - `gptca_confidence`
   - `gptca_status`
   - `gptca_rationale`
   - `gptca_ontology_id`
   - `gptca_proposed_label` (original suggestion when validation overrides the call)
4. Returns the updated `AnnData` and a structured report identical to the REST API payload.

### Configuration pointers
- Control the number of markers with `top_n_markers=` (defaults to 5).
- Set `result_prefix="custom"` to change column names (`custom_label`, etc.).
- Pass `marker_db_path="data/processed/marker_db.parquet"` or a pre-loaded DataFrame via `marker_db=` for custom knowledge bases.
- Supply an `Annotator` instance if you need specific settings (e.g., mock mode or custom API key).

## 3. Command-line workflow
Use the bundled CLI for non-interactive pipelines:
```bash
python -m gpt_cell_annotator.scanpy annotate data/demo/pbmc_demo.h5ad \
  --cluster-key leiden \
  --species "Homo sapiens" \
  --tissue Blood \
  --marker-db data/processed/marker_db.parquet \
  --output data/demo/pbmc_demo_annotated.h5ad \
  --summary-csv docs/reports/latest/pbmc_summary.csv \
  --summary-json docs/reports/latest/pbmc_report.json
```

Key flags:
- `--skip-recompute-markers` skips `scanpy.tl.rank_genes_groups` if you already computed rankings.
- `--prefix custom` writes `custom_label`, etc.
- `--summary-*` writes per-cluster CSV and the full JSON report (matches REST API output).

## 4. Troubleshooting
- **Missing marker DB**: run `python scripts/build_marker_db.py` or point `--marker-db` to a valid parquet file.
- **`rank_genes_groups` missing**: ensure Scanpy is installed; the helper will compute rankings unless `--skip-recompute-markers` is set.
- **Mock mode**: unset `OPENAI_API_KEY` to fall back to the heuristic annotator. Outputs remain deterministic for demos.
- **AnnData not updated**: confirm the `cluster_key` exists in `adata.obs` and that clusters have at least one marker gene after ranking.

## 5. What’s next?
- Extend notebooks by writing `report` JSON alongside `adata` to track evidence.
- Feed the annotations back into downstream Scanpy plots (e.g., `sc.pl.umap(adata, color="gptca_label")`).
- Combine with the benchmark toolkit to compare GPT Cell Annotator against reference-based approaches on your datasets.

Happy annotating!
