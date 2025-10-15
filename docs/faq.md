# Frequently Asked Questions

Looking for a quick answer? Start here and cross-reference the deeper guides linked throughout.

## Setup & Access

### Do I need an OpenAI API key?

Yes for live LLM annotations. Set `OPENAI_API_KEY` in your environment or `.env`. If the key is missing, expired, or rate-limited, the system automatically falls back to the heuristic mock annotator. You can force mock mode with `gca annotate ... --offline` or `Annotator(force_mock=True)`.

See: [`docs/install.md`](docs/install.md#command-line-interface).

### How do I run completely offline?

Use the bundled assets and mock annotator:

```bash
gca build-db --offline
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json
```

The first command copies packaged data into `~/.cache/gpt-cell-annotator`. Override with `GPT_CELL_ANNOTATOR_HOME`. Offline Docker usage is covered in [`docs/install.md`](docs/install.md#docker--compose).

### Where does the marker database live?

By default, under `~/.cache/gpt-cell-annotator/data/processed/`. Override with:

- `GPT_CELL_ANNOTATOR_DATA_DIR=/path/to/dir`
- `gca build-db --output-dir ./my_db`
- `annotate_anndata(..., marker_db_path="path/to/marker_db.parquet")`

## Annotation Behaviour

### Why do some clusters return “Unknown or Novel”?

The validation layer requires a minimum marker overlap (`VALIDATION_MIN_MARKER_OVERLAP`, default `2`). If the overlap is below the threshold or conflicting evidence is detected, the primary label is downgraded to `Unknown or Novel`. The original LLM suggestion is stored as `*_proposed_label` and `annotation.proposed_label`.

See: [`docs/operations.md`](docs/operations.md#validation-guardrails).

### How are gene synonyms and species differences handled?

Markers are normalised using `config/gene_synonyms.json` and the ortholog map in `data/orthologs/`. Provide the source species (`--species` or `dataset_context["species"]`) and the engine maps genes to the primary species (human by default). Update the synonym file or set `SYNONYM_ENABLE_ORTHOLOGS=false` for debugging.

See: [`docs/scanpy_integration.md`](docs/scanpy_integration.md#cross-species-workflows).

### Can I tune retrieval-augmented prompting (RAG)?

Yes. Environment variables:

- `RAG_ENABLED=false` disables database lookups.
- `RAG_TOP_K=<int>` and `RAG_MIN_OVERLAP=<int>` control candidate density.

Adjust these before running the API, CLI, or notebooks. The changes are picked up by `get_settings()`.

## Troubleshooting

### “Marker database not found” errors

Run `gca build-db` (or `python scripts/build_marker_db.py`). If you’re on a read-only system, copy `gpt_cell_annotator/_assets/data/processed/marker_db.parquet` to a writeable directory and set `GPT_CELL_ANNOTATOR_DATA_DIR`.

### CLI reports “rank_genes_groups missing”

Install the Scanpy extra (`pip install "gpt-cell-annotator[scanpy]"`) and avoid `--skip-recompute-markers` for AnnData files lacking differential-expression results. The helper will compute rankings automatically when Scanpy is available.

### API returns 500 with validation errors

Enable DEBUG logs (`LOG_LEVEL=DEBUG`) and inspect structured output. Common causes:

- Missing ontology IDs in LLM responses → ensure RAG is enabled.
- Malformed payloads → validate against `schemas/annotation_payload.schema.json`.

### Cached assets get stale

Delete the cache directory (`rm -rf ~/.cache/gpt-cell-annotator`) or pass `--output-dir` to `gca build-db` to regenerate assets in a fresh location. Remember to update `GPT_CELL_ANNOTATOR_HOME` to point to the new cache.

## Resources

- Installation & offline workflows: [`docs/install.md`](docs/install.md)
- Validation deep dive: [`docs/operations.md`](docs/operations.md)
- Benchmarks & evaluation: [`docs/benchmarks.md`](docs/benchmarks.md)
- Demo playbook: [`docs/demo.md`](docs/demo.md)

Still stuck? Open an issue with the failing command, environment details, and relevant logs.
