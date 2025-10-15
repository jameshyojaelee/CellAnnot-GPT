# Benchmark Datasets

This directory stores lightweight metadata for benchmark runs. Each JSON file
contains the dataset name, licensing information, optional download
instructions, and the list of clusters used for evaluation. The actual raw data
are **not** committed to the repository; follow the “source” links to obtain the
original atlases.

## Included samples

| File | Description | Source / License |
| --- | --- | --- |
| `pbmc_small.json` | Peripheral blood mononuclear cell toy subset (2 clusters). | Derived from 10x Genomics PBMC 3k, CC BY 4.0. |
| `brain_cortex_small.json` | Miniature human cortex subset (inhibitory vs excitatory neurons). | Based on Darmanis et al. 2015, CC BY 4.0. |
| `lung_epithelial_small.json` | Airway epithelial demo (AT2 vs Club cells). | Adapted from Travaglini et al. 2020, CC BY 4.0. |

To add a new dataset:

1. Create a `<name>.json` file mirroring the existing schema.
2. Include `dataset_context` with at least `species` and (optionally) `tissue`.
3. List clusters with `cluster_id`, `ground_truth`, and marker genes (top 5–10).
4. Document the licence and citation in this README or within the JSON file.
5. Store preprocessing notes or notebooks under `evaluation/datasets/<name>/`.

When running benchmarks, the CLI will look for all `*.json` files in this
directory unless a subset is specified explicitly.
