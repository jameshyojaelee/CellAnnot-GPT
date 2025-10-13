# CellAnnot-GPT

CellAnnot-GPT is an AI assistant that annotates single-cell RNA-seq clusters by combining curated marker knowledge with a large language model. It delivers evidence-backed cell type suggestions, confidence scoring, and validation tooling that reduces manual labeling time from days to minutes.

## Quick Start

```bash
# Clone and enter repo
git clone https://github.com/your-org/CellAnnot-GPT.git
cd CellAnnot-GPT

# Set up environment (Poetry example)
poetry install
poetry run poe build-marker-db   # or python scripts/build_marker_db.py

# Run FastAPI backend
poetry run uvicorn backend.api.main:app --reload

# Run Streamlit UI in another terminal
poetry run streamlit run frontend/streamlit_app.py
```

## Key Features

- Marker knowledge ingestion from PanglaoDB, CellMarker, and curated literature.
- Prompt-engineered LLM annotation engine with batch support and uncertainty handling.
- Validation layer cross-checking LLM outputs against marker databases.
- FastAPI REST endpoints plus Streamlit UI for interactive exploration.
- Benchmark toolkit for evaluating accuracy and tracking improvements.
- Deployment-ready configuration (Docker, CI/CD) and demo assets for presentations.

## Repository Layout

See directory tree in project root for component descriptions. Detailed architecture notes live in `docs/architecture.md`.
