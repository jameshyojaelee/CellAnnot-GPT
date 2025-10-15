# Benchmark Suite

GPT Cell Annotator ships with a lightweight benchmark harness that evaluates the
annotator against curated marker datasets and heuristic baselines. The suite is
designed to run offline (mock mode) yet can leverage live LLM calls when API
credentials are available.

## Datasets

Metadata lives under `evaluation/datasets/`:

| Dataset | Description | Notes |
| --- | --- | --- |
| `pbmc_small` | Peripheral blood mononuclear cell toy subset (B vs T cells). | Derived from 10x Genomics PBMC 3k (CC BY 4.0). |
| `brain_cortex_small` | Human cortex mini atlas (inhibitory/excitatory neurons + astrocytes). | Based on Darmanis et al. 2015 (CC BY 4.0). |
| `lung_epithelial_small` | Airway epithelial trio (AT2, Club, Basal). | Adapted from Travaglini et al. 2020 (CC BY 4.0). |

Each JSON file records the dataset context, citation, and a handful of clusters
with canonical marker genes and ground-truth labels. Extend the suite by adding
new JSON files and, optionally, preprocessing notebooks in
`evaluation/datasets/<name>/`.

## Running benchmarks

```bash
# Ensure marker DB is available
python scripts/build_marker_db.py

# Run all datasets in mock mode
python scripts/run_benchmarks.py --mock

# Restrict to specific datasets and output dir
python scripts/run_benchmarks.py \
  --datasets pbmc_small,brain_cortex_small \
  --output docs/reports \
  --baselines marker_overlap

# Convenience task via Poe
poetry run poe benchmarks
```

Outputs are written to `docs/reports/<timestamp>/` with the following structure:

- `<dataset>/<model>/metrics.json` – metrics per model (accuracy, macro F1, top-3
  recall, unknown rate, flag precision, time per cluster).
- `<dataset>/<model>/predictions.json` – per-cluster predictions with status and
  confidence.
- `<dataset>/<model>/*.png` – confusion matrix, precision/recall bars, and
  calibration plots.
- `<dataset>/README.md` – Markdown summary comparing models on that dataset.
- `summary.json` – aggregate view across datasets.
- `history.csv` – tabular archive suitable for dashboards.

By default, the baseline is a marker-overlap heuristic derived from the marker
database. Additional baselines can be registered by extending
`evaluation/benchmark_runner.py` (e.g., wrapping CellTypist predictions once
expression matrices are available).

## Interpreting metrics

- **Accuracy / Macro F1** – primary fidelity scores.
- **Top-3 recall** – fraction of clusters where the ground truth appears within
  the primary label plus the top two alternatives.
- **Unknown rate** – share of clusters flagged as `Unknown or Novel`.
- **Flag precision** – of flagged clusters, proportion where the predicted label
  disagreed with ground truth (i.e., helpful warnings).
- **Time per cluster** – average runtime, useful when comparing mock vs live LLM
  runs.

## Extending the suite

1. Add a dataset JSON file with `dataset_name`, `dataset_context`, `clusters`,
   and citation/licence fields.
2. Implement or register new baselines in `evaluation/benchmark_runner.py`.
3. Update this document with dataset descriptions.
4. Optionally wire the benchmark into CI or nightly automation via the `poe`
   task (see `pyproject.toml`).

Happy benchmarking!
