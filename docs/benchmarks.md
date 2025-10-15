# Benchmark Suite

GPT Cell Annotator ships with a lightweight benchmark harness that evaluates annotation quality against curated marker datasets and heuristic baselines. Everything runs offline by default; provide an API key to compare mock vs live LLM performance.

## Quickstart

| Step | Command | Purpose |
| --- | --- | --- |
| 1. Build assets | `gca build-db --offline` | Ensures marker DB is present for retrieval + validation. |
| 2. Run the full suite | `python scripts/run_benchmarks.py --mock` | Executes every dataset using the heuristic annotator. |
| 3. Compare live vs mock | `OPENAI_API_KEY=... python scripts/run_benchmarks.py --datasets pbmc_small --models live,mock` | Benchmarks both modes on PBMC. |
| 4. Inspect reports | `ls docs/reports/<timestamp>` | Review `summary.json` + per-dataset READMEs from the latest run. |

Additional convenience commands:

```bash
# Poetry shortcut defined in pyproject.toml
poetry run poe benchmarks

# Store reports under a custom directory and restrict datasets/baselines
python scripts/run_benchmarks.py \
  --datasets pbmc_small,brain_cortex_small \
  --baselines marker_overlap \
  --output docs/reports/$(date +%Y%m%d-%H%M)
```

## Datasets

Metadata lives under [`evaluation/datasets/`](../evaluation/datasets/):

| Dataset | Description | Notes |
| --- | --- | --- |
| `pbmc_small` | Peripheral blood mononuclear cell toy subset (B vs T cells). | Derived from 10x Genomics PBMC 3k (CC BY 4.0). |
| `brain_cortex_small` | Human cortex mini atlas (inhibitory/excitatory neurons + astrocytes). | Based on Darmanis et al. 2015 (CC BY 4.0). |
| `lung_epithelial_small` | Airway epithelial trio (AT2, Club, Basal). | Adapted from Travaglini et al. 2020 (CC BY 4.0). |
| `mouse_pbmc_small` | Mouse PBMC subset (B vs T cells; ortholog demo). | Adapted from 10x Genomics mouse PBMC (CC BY 4.0). |

Each JSON file records dataset context, citation, and clusters with canonical markers + ground truth labels. Extend the suite by adding new JSON files (and optional preprocessing notebooks) under `evaluation/datasets/<name>/`.

## Output Structure

Benchmark runs emit timestamped directories under `docs/reports/`:

- `<dataset>/<model>/metrics.json` – accuracy, macro F1, top-3 recall, unknown rate, flag precision, and runtime metrics.
- `<dataset>/<model>/predictions.json` – per-cluster outputs with rationale, confidence, and validation status.
- `<dataset>/<model>/*.png` – confusion matrices, calibration plots, and marker overlap charts (generated when matplotlib is available).
- `<dataset>/README.md` – human-readable summary comparing models on that dataset.
- `summary.json` – aggregate metrics across datasets.
- `history.csv` – append-only history suitable for dashboards.

Pair the reports with the UI or `docs/demo_assets/pbmc_annotation.json` for live walkthroughs.

## Interpreting Metrics

- **Accuracy / Macro F1** – ground-truth alignment across clusters.
- **Top-3 recall** – likelihood the true label appears in the primary label or the top two alternatives.
- **Unknown rate** – proportion of clusters flagged as `Unknown or Novel` after validation.
- **Flag precision** – fraction of flagged clusters where the prediction disagreed with ground truth (higher is better).
- **Time per cluster** – average runtime per cluster; compare live vs mock or alternative backends.

Use these figures to calibrate guardrails (see [`docs/operations.md`](docs/operations.md)) and to demonstrate gains over reference-based methods (see comparison in [`docs/why.md`](docs/why.md#state-of-the-field)).

## Extending & Customising

1. Add a dataset JSON file with `dataset_name`, `dataset_context`, `clusters`, and citation/licence fields.
2. Implement additional baselines in [`evaluation/benchmark_runner.py`](../evaluation/benchmark_runner.py) (e.g., wrap CellTypist predictions).
3. Register the baseline in `scripts/run_benchmarks.py` and update this doc’s dataset table.
4. Wire the benchmark task into CI or nightly jobs via the `poe` task or a GitHub Actions workflow.

Happy benchmarking!
