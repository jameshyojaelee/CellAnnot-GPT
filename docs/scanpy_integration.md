# Scanpy Integration Guide

Leverage GPT Cell Annotator inside Scanpy notebooks, pipelines, and headless batch jobs. The integration ships with both a rich Python API (`annotate_anndata`) and a CLI (`gca scanpy`) so you can stay in your preferred workflow.

## Quickstart Checklist

| Step | Command | Notes |
| --- | --- | --- |
| 1. Install extras | `pip install "gpt-cell-annotator[scanpy]"` | Adds Scanpy + AnnData dependencies. |
| 2. Materialise assets | `gca build-db --offline` | Uses bundled demo sources; see [`scripts/build_marker_db.py`](../scripts/build_marker_db.py) for custom runs. |
| 3. Annotate markers from a CSV | `gca annotate data/demo/pbmc_markers.csv --species "Homo sapiens" --offline` | Produces JSON/CSV reports you can merge back into AnnData. |
| 4. Port results into AnnData | Follow the code below or the walkthrough in [`docs/demo.md`](docs/demo.md) | Demonstrates notebook + UI steps end-to-end. |

> Need end-to-end setup guidance? Start with `docs/getting_started.md` and the CLI primer in `docs/install.md`.

## Annotate Within a Notebook

```python
import scanpy as sc
from gpt_cell_annotator import annotate_anndata

adata = sc.read_h5ad("data/demo/pbmc_demo.h5ad")

# Ensure differential expression is computed (will be done automatically if missing)
sc.tl.rank_genes_groups(adata, groupby="leiden", n_genes=5, method="wilcoxon")

adata, report = annotate_anndata(
    adata,
    cluster_key="leiden",
    species="Homo sapiens",
    tissue="Peripheral blood",
    top_n_markers=5,
)

adata.obs[["gptca_label", "gptca_confidence", "gptca_status"]].head()
```

### What the helper does

1. Computes `rank_genes_groups` if the AnnData object lacks marker rankings.
2. Normalises gene symbols and builds the batch payload.
3. Calls the annotation engine (mock or live depending on `OPENAI_API_KEY`).
4. Cross-checks results against the marker DB and writes columns such as `gptca_label`, `gptca_status`, `gptca_rationale`, `gptca_canonical_markers`, and `gptca_mapping_notes`.
5. Returns the updated `AnnData` plus a structured report identical to the REST API schema (`docs/api_reference.md`).

### Common Parameters

- `top_n_markers`: choose how many markers per cluster are sent to the LLM (default `5`).
- `result_prefix`: customise output column prefixes (e.g., `"ca2025_label"`).
- `marker_db_path` / `marker_db`: point to custom knowledge bases (Parquet DataFrame).
- `annotator`: inject your own `Annotator` (e.g., forcing `force_mock=True` for validation runs).
- `compute_rankings=False`: skip automatic `scanpy.tl.rank_genes_groups` if you already have rankings.

## Command-Line Workflow

The CLI mirrors the notebook flow and is ideal for automation, Nextflow/Snakemake steps, or CI smoke tests.

```bash
gca scanpy annotate data/demo/pbmc_demo.h5ad \
  --cluster-key leiden \
  --species "Homo sapiens" \
  --tissue "Peripheral blood" \
  --summary-json reports/pbmc_report.json \
  --summary-csv reports/pbmc_summary.csv \
  --offline
```

Key flags:

- `--marker-db` – override the Parquet file. Defaults to the cached `~/.cache/gpt-cell-annotator/data/processed/marker_db.parquet`.
- `--skip-recompute-markers` – prevent recomputing `rank_genes_groups`.
- `--prefix custom` – rename output columns (`custom_label`, etc.).
- `--offline` / `--mock` – force heuristic annotations and avoid external API calls.

For reference implementations see [`docs/demo.md`](docs/demo.md) for the live demo script and [`scripts/run_benchmarks.py`](../scripts/run_benchmarks.py) for a headless example that feeds annotation outputs into evaluation metrics.

## Validation & Guardrails

Results returned by `annotate_anndata` or `gca scanpy` run through the same validation layer described in `docs/operations.md#validation-guardrails`. Expect:

- Automatic downgrade to `Unknown or Novel` when marker overlap falls below `VALIDATION_MIN_MARKER_OVERLAP`.
- `*_mapping_notes` columns containing ortholog conversions and missed genes.
- Structured warnings persisted in the JSON report and `*_status` fields (e.g., `flagged_low_support`, `flagged_species_mismatch`).

Tune thresholds via `.env` or command-line environment variables before invoking the CLI:

```bash
export VALIDATION_MIN_MARKER_OVERLAP=3
export CONFIDENCE_OVERLAP_HIGH=4
gca scanpy annotate ...
```

## Cross-Species Workflows

- Provide the appropriate `--species` or `dataset_context["species"]` (e.g., `"Mus musculus"`).
- Ortholog expansion uses files in `data/orthologs/`; customise via `ORTHOLOG_MAPPING_PATH`.
- Inspect `gptca_canonical_markers` to see the humanised markers used during validation.
- Disable ortholog handling with `SYNONYM_ENABLE_ORTHOLOGS=false` for debugging or when working with custom gene IDs.

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| `FileNotFoundError: marker_db.parquet` | Run `gca build-db` or specify `--marker-db`. |
| `ValueError: rank_genes_groups` missing | Install the Scanpy extra and avoid `--skip-recompute-markers`. |
| Mock annotations when API key is set | Confirm `OPENAI_API_KEY` is exported in the environment invoking the CLI/notebook. |
| No markers written to `*_canonical_markers` | Ensure clusters contain at least one ranked gene above `top_n_markers`. |

## Further Reading

- `docs/install.md` – package installation and offline usage.
- `docs/getting_started.md` – full environment tutorial including Scanpy workflows.
- `docs/operations.md` – guardrails, validation thresholds, and observability tips.
- `docs/benchmarks.md` – comparing Scanpy-integrated runs with baseline models.

Happy annotating!
