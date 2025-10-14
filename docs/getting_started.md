# Getting Started with CellAnnot-GPT

Welcome! This guide walks you through setting up the development environment, running tests, and trying the core workflows.

## 1. Prerequisites
- Python 3.11
- Poetry (recommended) or pip
- Optional: Docker for containerized workflows

## 2. Environment Setup
```bash
# Clone the repo
git clone https://github.com/your-org/CellAnnot-GPT.git
cd CellAnnot-GPT

# Install dependencies with Poetry
poetry install

# Alternatively, using pip and venv
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Copy environment template and fill values
cp .env.template .env
```

## 3. Run the stack locally
```bash
# Build marker knowledge base (downloads remote sources by default)
poetry run python scripts/build_marker_db.py

# Offline mode → stick to bundled demo assets
poetry run python scripts/build_marker_db.py --local-only

# Verify downloaded snapshots against checksums in config/marker_sources.yaml
poetry run python scripts/build_marker_db.py --verify-checksums

# Start FastAPI backend
poetry run uvicorn backend.api.main:app --reload

# In another terminal, launch the Streamlit UI
poetry run streamlit run frontend/streamlit_app.py
```
Visit `http://localhost:8501` to explore the dashboard. Each entry in `config/marker_sources.yaml` records a live download URL, checksum, and version tag. Pass `--local-only` to reuse the bundled demo assets or `--verify-checksums` to enforce snapshot pinning during ingestion.

## 4. Run tests
```bash
poetry run pytest            # full test suite
poetry run ruff check        # lint
poetry run mypy backend      # type checks
```

## 5. Sample workflows
- **Batch annotation**: upload `data/demo/pbmc_markers.csv` via the UI and review results.
- **Single cluster**: switch to the “Single Cluster” tab, enter markers (e.g., `MS4A1, CD79A`) to get instant predictions.
- **Benchmarks**: `poetry run python scripts/run_benchmarks.py --mock` will generate reports under `docs/reports/<date>/`.

## 6. Useful documentation
- Architecture: see `docs/architecture.md`
- Deployment: see `docs/deployment.md`
- Operations & logging: see `docs/operations.md`
- Demo playbook: see `docs/demo.md`

## 7. Troubleshooting
- Missing markers DB → rerun the build script.
- Mock mode warning → set `OPENAI_API_KEY` in `.env` to enable live annotations.
- Redis cache disabled → set `REDIS_URL` if you need caching.

Happy hacking!
