# GPT Cell Annotator

GPT Cell Annotator is an AI assistant that annotates single-cell RNA-seq clusters by combining curated marker knowledge with a large language model. It delivers evidence-backed cell type suggestions, confidence scoring, and validation guardrails so teams can move from raw clusters to trusted labels in minutes.

## Why It Matters

- Manual scRNA-seq annotation is slow and inconsistent; experts triage clusters by cross-referencing markers, literature, and ontologies — hours to days per dataset.
- Existing tools (reference-mapping like SingleR/CellTypist, latent-transfer like scANVI) work best with high-quality reference atlases; they struggle on novel cell states, rare types, or out-of-distribution data.
- An LLM assistant adds flexible knowledge integration and explanations, but must be guarded against hallucinations via validation and “unknown” handling.

See `docs/why.md` for background and motivation.

## Quick Start

### 1. Install

```bash
pip install "gpt-cell-annotator[api]"
# Optional extras:
#   [ui]     -> Streamlit dashboard + plotting stack
#   [scanpy] -> Scanpy integration helpers
```

### 2. Annotate markers from the CLI

```bash
# Fully offline demo (bundled assets + heuristic mock annotator)
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json

# Rebuild the marker database locally
gca build-db --offline --output-dir ~/.cache/gca/db
```

### 3. Serve the API or UI

```bash
# API (use --offline to disable OpenAI calls)
OPENAI_API_KEY=sk-... gca api --host 0.0.0.0 --port 8000

# Streamlit UI (requires [ui] extra)
streamlit run frontend/streamlit_app.py
```

Sample API calls (after the API is running on `localhost:8000`):

```bash
curl -s http://127.0.0.1:8000/health

curl -s -X POST http://127.0.0.1:8000/annotate_cluster \
  -H "Content-Type: application/json" \
  -d '{"cluster": {"cluster_id": "0", "markers": ["MS4A1", "CD79A"]},
       "dataset_context": {"species": "Homo sapiens", "tissue": "Blood"}}'
```

### 4. Containers

```bash
# Build the multi-stage image and launch the offline stack
docker compose up --build

# One-off annotation using the CLI entrypoint
docker run --rm -v $PWD/data:/data gpt-cell-annotator annotate data/demo/pbmc_markers.csv --offline
```

> Offline mode relies on the bundled heuristic annotator and local marker database copies. Live LLM annotations require `OPENAI_API_KEY`.

See `docs/install.md` for detailed installation scenarios.

## Quick Demo

```bash
gca build-db --offline
gca annotate data/demo/pbmc_markers.csv --offline --species "Homo sapiens"
```

Expected console output:

```
Cluster  Primary Label     Confidence  Status     Warnings
0        B cell            High        supported  -
5        Unknown or Novel  Low         flagged    low_marker_overlap
```

Add `--out-json demo_annotations.json` to store the structured report (same schema as the REST API). Continue with the UI or Scanpy tutorials in [`docs/getting_started.md`](docs/getting_started.md#guided-tutorials).

## Key Features

[![Getting Started](https://img.shields.io/badge/docs-getting_started-blue)](docs/getting_started.md)
[![Install Guide](https://img.shields.io/badge/docs-install-blueviolet)](docs/install.md)
[![API Reference](https://img.shields.io/badge/docs-api_reference-green)](docs/api_reference.md)
[![Scanpy Integration](https://img.shields.io/badge/docs-scanpy_integration-teal)](docs/scanpy_integration.md)
[![Operations](https://img.shields.io/badge/docs-operations-purple)](docs/operations.md)
[![Benchmarks](https://img.shields.io/badge/docs-benchmarks-red)](docs/benchmarks.md)
[![FAQ](https://img.shields.io/badge/docs-faq-lightgrey)](docs/faq.md)
[![Roadmap](https://img.shields.io/badge/docs-roadmap-yellow)](docs/roadmap.md)

**Workflow cheat sheet**
- Build the marker knowledge base from the sources listed in `config/marker_sources.yaml` with `gca build-db`; outputs land in `data/processed/`.
- Start the FastAPI backend (`gca api`) so it can load the marker DB and expose `/annotate_cluster` / `/annotate_batch`.
- Upload cluster markers (e.g., `data/demo/pbmc_markers.csv`) via the Streamlit UI or call the API to receive JSON annotations you can store alongside Scanpy results.
- Drop into notebooks with `annotate_anndata` or run `gca scanpy annotate` for batch pipelines.
- Tune guardrails with `VALIDATION_MIN_MARKER_OVERLAP`, `VALIDATION_FORCE_UNKNOWN_ON_FAIL`, `CONFIDENCE_OVERLAP_MEDIUM`, and `CONFIDENCE_OVERLAP_HIGH` to balance conservatism vs. recall.
- Cross-species? Supply `species` (e.g., `Mus musculus`)—ortholog tables in `data/orthologs/` map markers back to the human-centric knowledge base and surface mapping notes in the UI/API.

- Marker knowledge ingestion from PanglaoDB, CellMarker, and curated literature.
- Prompt-engineered LLM annotation engine with batch support and uncertainty handling.
- Validation layer cross-checking LLM outputs against marker databases.
- FastAPI REST endpoints plus a Streamlit UI for interactive exploration, run-to-run comparison, and knowledge-linked warnings.
- Benchmark toolkit for evaluating accuracy and tracking improvements.
- Deployment-ready configuration (Docker, CI/CD) and demo assets for presentations.

## Repository Layout

See directory tree in project root for component descriptions. Detailed architecture notes live in `docs/architecture.md`.
