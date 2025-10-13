# Benchmarking CellAnnot-GPT

This guide explains how to evaluate the annotator against curated datasets.

## Datasets
- Store benchmark configs under `evaluation/datasets/`. Each JSON file should contain:
  ```json
  {
    "dataset_name": "pbmc_small",
    "dataset_context": {"species": "Homo sapiens", "tissue": "Peripheral blood"},
    "clusters": [
      {"cluster_id": "0", "ground_truth": "B cell", "markers": ["MS4A1", "CD79A"]}
    ]
  }
  ```
- Include `ground_truth` labels for scoring.

## Running benchmarks
```bash
poetry run python scripts/run_benchmarks.py --datasets-dir evaluation/datasets --output-dir docs/reports
```

Options:
- `--datasets` specific files to run (space-separated)
- `--mock` forces mock annotator (no OpenAI key required)

Outputs are written to `docs/reports/YYYYMMDD/<dataset>.{json,md}`.

## Interpreting reports
- JSON: structured scores (accuracy, macro F1, per-class metrics, predictions).
- Markdown: formatted summary plus text confusion matrix for quick review.

## CI automation
- Workflow `ci.yml` includes a scheduled job (`cron`) to run benchmarks nightly.
- Reports are uploaded as GitHub Actions artifacts for auditing.
- Ensure secrets (`OPENAI_API_KEY`) are configured in repository settings if running live mode.

