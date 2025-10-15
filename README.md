# GPT Cell Annotator

GPT Cell Annotator is an AI assistant that annotates single-cell RNA-seq clusters by combining curated marker knowledge with a large language model. It delivers evidence-backed cell type suggestions, confidence scoring, and validation tooling that reduces manual labeling time from days to minutes.

## Why It Matters

- Manual scRNA-seq annotation is slow and inconsistent; experts triage clusters by cross-referencing markers, literature, and ontologies — hours to days per dataset.
- Existing tools (reference-mapping like SingleR/CellTypist, latent-transfer like scANVI) work best with high-quality reference atlases; they struggle on novel cell states, rare types, or out-of-distribution data.
- An LLM assistant adds flexible knowledge integration and explanations, but must be guarded against hallucinations via validation and “unknown” handling.

See `docs/why.md` for background and motivation.

## Quick Start

```bash
# Clone and enter repo
git clone https://github.com/your-org/GPT-cell-annotator.git
cd GPT-cell-annotator

# Set up environment (Poetry example)
poetry install
poetry run poe build-marker-db   # or python scripts/build_marker_db.py

# Run FastAPI backend
poetry run uvicorn backend.api.main:app --reload

# Run Streamlit UI in another terminal
poetry run streamlit run frontend/streamlit_app.py

# Docker (optional)
docker compose up --build

# Sample API calls (after backend is running on localhost:8000)
curl -s http://127.0.0.1:8000/health

curl -s -X POST http://127.0.0.1:8000/annotate_cluster \
  -H "Content-Type: application/json" \
  -d '{
        "cluster": {"cluster_id": "0", "markers": ["MS4A1", "CD79A"]},
        "dataset_context": {"species": "Homo sapiens", "tissue": "Blood"}
      }'

curl -s -X POST http://127.0.0.1:8000/annotate_batch \
  -H "Content-Type: application/json" \
  -d '{
        "clusters": [
          {"cluster_id": "0", "markers": ["MS4A1"]},
          {"cluster_id": "1", "markers": ["CD3E"]}
        ],
        "dataset_context": {"species": "Homo sapiens"}
      }'
```

## Key Features

[![Getting Started](https://img.shields.io/badge/docs-getting_started-blue)](docs/getting_started.md)
[![API Reference](https://img.shields.io/badge/docs-api_reference-green)](docs/api_reference.md)
[![Operations](https://img.shields.io/badge/docs-operations-purple)](docs/operations.md)
[![Workflow Overview](https://img.shields.io/badge/docs-workflow-orange)](docs/getting_started.md#workflow-at-a-glance)
[![Scanpy Integration](https://img.shields.io/badge/docs-scanpy_integration-teal)](docs/scanpy_integration.md)

**Workflow cheat sheet**
- Build the marker knowledge base from the sources listed in `config/marker_sources.yaml`; outputs land in `data/processed/`.
- Start the FastAPI backend so it can load the marker DB and expose `/annotate_cluster` / `/annotate_batch`.
- Upload cluster markers (e.g., `data/demo/pbmc_markers.csv`) via the Streamlit UI or call the API to receive JSON annotations you can store alongside Scanpy results.
- Drop into notebooks with `annotate_anndata` or run `python -m gpt_cell_annotator.scanpy annotate` for batch pipelines.

- Marker knowledge ingestion from PanglaoDB, CellMarker, and curated literature.
- Prompt-engineered LLM annotation engine with batch support and uncertainty handling.
- Validation layer cross-checking LLM outputs against marker databases.
- FastAPI REST endpoints plus a Streamlit UI for interactive exploration, run-to-run comparison, and knowledge-linked warnings.
- Benchmark toolkit for evaluating accuracy and tracking improvements.
- Deployment-ready configuration (Docker, CI/CD) and demo assets for presentations.

## Repository Layout

See directory tree in project root for component descriptions. Detailed architecture notes live in `docs/architecture.md`.
