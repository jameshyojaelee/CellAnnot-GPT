# Installation Guide

This guide walks through installing GPT Cell Annotator from Python packages, containers, and source. It also covers offline usage and asset management.

## Prerequisites

- Python 3.11
- `pip >= 23.0`
- ~2 GB of free disk space for optional extras (UI, Scanpy)

## Install with `pip`

```bash
# Base CLI + annotation engine
pip install gpt-cell-annotator

# Include FastAPI + Redis support for the HTTP API
pip install "gpt-cell-annotator[api]"

# Include Streamlit UI and plotting components
pip install "gpt-cell-annotator[ui]"

# Add Scanpy integration helpers (if you plan to use `gca scanpy`)
pip install "gpt-cell-annotator[scanpy]"

# Bring in development tooling (Ruff, mypy, pytest, twine)
pip install "gpt-cell-annotator[dev]"
```

The package bundles demo assets, marker databases, and default configuration. The first time you run a command the assets are copied into a cache directory:

- Default cache: `~/.cache/gpt-cell-annotator`
- Override: set `GPT_CELL_ANNOTATOR_HOME=/path/to/cache`
- Force a specific marker database: `GCA_MARKER_DB_PATH=/mnt/atlas/marker_db.parquet`
- Persist annotation responses offline: `GCA_CACHE_DIR=~/.cache/gca/annotations`

## Command Line Interface

The `gca` CLI is available after installation:

```bash
# Annotate demo markers completely offline (uses mock LLM backend)
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json

# Rebuild the marker database (local sources only)
gca build-db --offline --output-dir ~/.cache/gca/db

# Forward to the Scanpy helper CLI
gca scanpy annotate demo.h5ad --species "Homo sapiens" --cluster-key leiden \
  --batch-size 16 --concurrency 2 --cache-dir ~/.cache/gca/annotations

gca scanpy validate demo.h5ad --species "Homo sapiens" --cluster-key leiden \
  --label-column curated_label

# Launch the FastAPI server
gca api --offline --host 0.0.0.0 --port 8000
```

Passing `--offline` or `--mock` forces the heuristic annotator and prevents external HTTP requests. Without the switch, the CLI will use your configured `OPENAI_API_KEY` for live annotations.

## Install for Seurat (R)

- Ensure R ≥ 4.3, Pandoc ≥ 2.11, and development headers for libcurl/openssl/xml2 are installed. On macOS use the Xcode command line tools; on Windows install [Rtools42](https://cran.r-project.org/bin/windows/Rtools/).
- Install directly from GitHub with pak (preferred) or remotes:

```r
if (!require("pak")) install.packages("pak")
pak::pkg_install("github::jameshyojaelee/CellAnnot-GPT@main?subdir=clients/r/gptcellannotator")
# remotes alternative
remotes::install_github("jameshyojaelee/CellAnnot-GPT", subdir = "clients/r/gptcellannotator", build = TRUE)
```

- Use the helper scripts for reproducible builds:
  - `Rscript clients/r/scripts/install_github.R` installs from the default branch (override with `Rscript clients/r/scripts/install_github.R v0.1.0`).
  - `Rscript clients/r/scripts/install_local_tarball.R` builds and installs a source tarball while caching dependencies.

- After installation, surface the pkgdown reference at <https://jameshyojaelee.github.io/CellAnnot-GPT/r> and run `vignette("annotate-seurat", package = "gptcellannotator")` to confirm the Seurat workflow renders.

- Offline CLI fallback relies on the Python package; install it alongside R with `pip install gpt-cell-annotator` and ensure `gca` is discoverable via `Sys.which("gca")`.

## Docker & Compose

The repository ships with a multi-stage Dockerfile that produces slim runtime images and a Docker Compose configuration tuned for offline demos.

```bash
# Build the container image
docker build -t gpt-cell-annotator .

# Run a one-off annotation
docker run --rm -v $PWD/data:/data gpt-cell-annotator annotate data/demo/pbmc_markers.csv --offline

# Launch the API + Redis stack (offline by default)
docker compose up --build
```

Volumes:

- `/data` inside the container maps to bundled assets and processed databases.
- Override by setting `GPT_CELL_ANNOTATOR_HOME` or `GPT_CELL_ANNOTATOR_DATA_DIR`.

## Working from Source

```bash
git clone https://github.com/your-org/GPT-cell-annotator.git
cd GPT-cell-annotator

python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e ".[api,ui,scanpy]"

# Smoke test the CLI
gca annotate data/demo/pbmc_markers.csv --offline
```

## Offline Checklist

- Set `--offline` / `--mock` on CLI commands.
- Ensure `GPT_CELL_ANNOTATOR_HOME` points to a writeable location (defaults to `~/.cache/gpt-cell-annotator`).
- For Docker/Compose, mount a persistent volume to `/data` so cached assets survive container rebuilds.
- When running the API offline, expect heuristic annotations; results from the live LLM require an `OPENAI_API_KEY`.
