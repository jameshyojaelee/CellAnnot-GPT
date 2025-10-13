# CellAnnot-GPT Build Prompt Pack

Use this sequence of GPT-5 prompts to generate the full CellAnnot-GPT project. Feed them one at a time, review outputs, and copy resulting code/files into your workspace. Replace bracketed placeholders with your specifics before sending.

## P0 – Vision & Constraints Alignment
```text
You are assisting with the CellAnnot-GPT project, an AI assistant that annotates single-cell RNA-seq clusters.

Deliver:
1. A concise architecture overview covering backend services, data layer, LLM interface, and UI.
2. Key technical decisions (programming languages, frameworks, hosting targets) and the rationale.
3. A list of core files/folders we must implement in the 4-week MVP.
4. Risks + mitigations (accuracy, latency, privacy).

Project constraints:
- Primary language: Python 3.11 for backend; Streamlit or FastAPI + simple front end.
- LLM: start with OpenAI GPT-4o API, but allow swapping to local model later.
- Knowledge sources: PanglaoDB, CellMarker, curated literature snippets.

Respond with a bullet summary ready for copy/paste into the project README "Overview" section.
```

## P1 – Repository Skeleton & Documentation
```text
Act as a senior software engineer setting up the CellAnnot-GPT repo.

Tasks:
1. Output the directory tree with brief purpose notes for each folder/file.
2. Provide contents for:
   - README.md (project intro, quick start, feature list).
   - docs/architecture.md (diagram description + component responsibilities).
   - .gitignore suitable for Python, Streamlit, and VSCode.
3. Include instructions for creating these files via shell commands.

Return the directory tree first, then fenced code blocks for each file.
```

## P2 – Environment & Dependency Setup
```text
Goal: Produce configuration files to manage dependencies.

Deliverables:
1. pyproject.toml using Poetry with sections for main deps (openai, pandas, numpy, anndata, fastapi or streamlit, uvicorn, scikit-learn, rich) and dev deps (pytest, black, mypy, ruff).
2. requirements.txt equivalent (for users without Poetry).
3. setup scripts: makefile targets (install, format, lint, test) inside Makefile.
4. Optional: environment.yml for conda users.
Return each file in separated fenced blocks.
```

## P3 – Knowledge Base Data Pipeline
```text
Design the data ingestion layer for marker gene references.

Provide:
1. Source summary table (PanglaoDB, CellMarker, curated literature).
2. Python module `backend/data_ingest/marker_loader.py` that downloads/parses raw CSV/JSON into a normalized parquet/SQLite store.
3. Tests in `tests/test_marker_loader.py` using fixtures/mocks (no real network calls).
4. CLI entry point `scripts/build_marker_db.py` to run the pipeline.

Output code blocks for each file. Ensure docstrings explain data schema.
```

## P4 – LLM Annotation Engine
```text
Implement the core annotation engine.

Produce:
1. `backend/llm/prompts.py` with reusable prompt templates (system, single cluster, batch, uncertainty) as Python functions.
2. `backend/llm/annotator.py` containing Annotator class with methods:
   - `annotate_cluster(cluster_payload)`
   - `annotate_batch(list_of_clusters)`
   - internal helper that calls OpenAI API with retry + rate limiting.
3. Config file `config/settings.py` using pydantic BaseSettings for API keys, model name, temperature, etc.
4. Unit tests `tests/test_annotator.py` mocking OpenAI client.

Return code with comments for non-obvious logic.
```

## P5 – Reasoning Guardrails & Validation
```text
Extend the engine with validation layer.

Deliver:
1. `backend/validation/crosscheck.py` that compares LLM output against marker DB (ontology alignment, contradictory marker detection).
2. `backend/validation/report.py` generating structured JSON + human-readable summary of annotations, warnings, and unknown clusters.
3. Example JSON schema definition in `schemas/annotation_result.schema.json`.
4. Tests `tests/test_validation.py` covering pass/fail scenarios.
```

## P6 – API Layer (FastAPI Option)
```text
Build a REST API around the annotation engine.

Provide:
1. `backend/api/main.py` with FastAPI app exposing endpoints:
   - `POST /annotate_cluster`
   - `POST /annotate_batch`
   - `GET /health`
2. Pydantic request/response models in `backend/api/models.py`.
3. Middleware for request logging + timing in `backend/api/middleware.py`.
4. `tests/test_api.py` using FastAPI TestClient and mocked annotator.
5. Update README quick start with curl examples.
```

## P7 – Streamlit Front-End (Alternative or Companion UI)
```text
Create a simple Streamlit UI consuming the API or local annotator.

Deliver:
1. `frontend/streamlit_app.py` with pages: Upload markers, Batch annotate, Review results.
2. Utility module `frontend/utils.py` handling API calls, caching, status indicators.
3. Instructions for running locally (`streamlit run frontend/streamlit_app.py`).
4. UI-focused tests using streamlit-testing-tools or minimal snapshot approach (describe strategy if tooling limited).
```

## P8 – Benchmark & Evaluation Toolkit
```text
Implement benchmarking utilities.

Provide:
1. `evaluation/benchmark_runner.py` that loads ground-truth datasets, runs annotations, computes accuracy/F1.
2. `evaluation/report_templates.py` to render markdown or HTML reports summarizing metrics, failure cases, confusion matrix (text-based if plotting unavailable).
3. Sample config `evaluation/datasets/pbmc_small.json` referencing test files.
4. Tests `tests/test_benchmark_runner.py` using synthetic data.
```

## P9 – Deployment & Ops
```text
Prepare deployment assets.

Deliver:
1. Dockerfile (multi-stage: builder + runtime) with Poetry install.
2. docker-compose.yml wiring API + optional Postgres/Redis if needed.
3. GitHub Actions workflow `.github/workflows/ci.yml` running lint, tests, and building Docker image.
4. docs/deployment.md explaining environments (dev/stage/prod) and secrets management.
```

## P10 – Admissions Demo Assets
```text
Generate materials for the 3-minute admissions demo.

Provide:
1. Updated script tailored to the current build.
2. Slide outline (titles + bullet points) for 3 slides: Problem, Live Demo, Impact.
3. Demo checklist (datasets, commands, failure fallback plan).

Deliver in Markdown sections ready to copy into docs/demo.md.
```

---

### Workflow Tips
- After each GPT-5 response, implement the files exactly, run lint/tests, and adjust prompts if gaps remain.
- Track decisions in docs/architecture.md so future prompts stay aligned.
- When swapping GPT models, reuse the same prompt structures to keep consistency.
