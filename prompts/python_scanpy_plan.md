# Mega Prompt: Elevate GPT Cell Annotator Python/Scanpy Experience

You are GPT-5 Codex acting as the Python lead. Your mission is to make GPT Cell Annotator the best-in-class solution for Scanpy/AnnData users and to polish packaging for broad distribution.

## Goals
- Harden and document the existing `annotate_anndata` API and `gca scanpy` CLI so they feel native to Scanpy workflows.
- Publish a polished Python distribution on PyPI with extras (`[scanpy]`, `[ui]`, `[api]`) and ensure consistent Docker + CLI experiences.
- Provide notebooks, benchmarks, and automation that showcase superiority over GPTCelltype.
- Add telemetry, caching, and offline ergonomics that reinforce reliability.

## Constraints / Assumptions
- Source tree: `gpt_cell_annotator` Python package (pyproject/Poetry managed).
- Current integration located in `gpt_cell_annotator/scanpy.py` with tests in `tests/test_scanpy_integration.py`.
- Must maintain compatibility with AnnData ≥ 0.9, Scanpy ≥ 1.9.
- Assume backend services and marker DB assets follow existing schema (`schemas/annotation.py`).
- CI uses GitHub Actions; prefer `uv` or `poetry run` for commands.

## Deliverables
1. Packaging & Distribution
   - Update `pyproject.toml` to define extras (`scanpy`, `ui`, `dev`), entry points, and metadata.
   - Automated build pipeline producing sdist + wheels, publishable via `uv build` or `poetry`.
   - Release checklist in `RELEASING.md` covering PyPI upload, Docker tagging, checksum verification.
2. Scanpy API Enhancements
   - Expand `annotate_anndata` to support asynchronous batching, cached marker DB loading, configurable guardrail thresholds, and progress logging.
   - Add convenience wrappers: `annotate_rank_genes`, `annotate_from_markers`.
   - Provide typed return objects (e.g., `DatasetReportModel`) and ensure mypy coverage.
3. CLI Improvements
   - Extend `gca scanpy annotate` to accept Loom/h5ad, chunked processing, species presets, offline toggles, and JSON report output path.
   - Add `gca scanpy validate` command for guardrail-only checks.
4. Documentation & Tutorials
   - Refresh `docs/scanpy_integration.md` with step-by-step guides, offline tips, and troubleshooting.
   - Create Jupyter notebooks under `notebooks/` demonstrating PBMC workflow, benchmarking vs. GPTCelltype, and guardrail visualization.
   - Update README badges/instructions (`pip install "gpt-cell-annotator[scanpy]"`).
5. Testing & QA
   - Expand unit tests (`tests/test_scanpy_integration.py`) to cover new features, caching, CLI options, and validation logic.
   - Integration tests with synthetic AnnData fixtures (store under `tests/fixtures/`).
   - Add end-to-end smoke test in CI running CLI with offline heuristic.
6. Telemetry & Operations
   - Emit structured logs (using `structlog`) with request IDs.
   - Optional counters/metrics hooks (Prometheus) behind extras.
   - Document configuration via env vars (`GCA_MARKER_DB_PATH`, `GCA_CACHE_DIR`).
7. Developer Tooling
   - Ensure Ruff/mypy/pytest all run via `make test` and CI.
   - Add pre-commit hooks for formatting and static analysis.

## Execution Steps
1. **Discovery & Inventory**
   - Review `gpt_cell_annotator/scanpy.py` and existing tests to understand payload flow.
   - Audit packaging metadata and CLI entry points in `pyproject.toml`.
   - List current docs and identify gaps vs. new goals.
2. **Packaging Upgrades**
   - Modify `pyproject.toml` extras and dependencies; ensure optional dependencies stay optional.
   - Configure build backend (e.g., `poetry-core` or `setuptools.build_meta`) for wheel generation.
   - Implement `make release` script (build, twine upload dry-run).
3. **API Enhancements**
   - Refactor `annotate_anndata` for batch streaming, caching, and guardrail overrides.
   - Add docstrings, type hints, and ensure compatibility with dataclasses/pydantic models.
   - Update `__all__` exports and module docs.
4. **CLI Expansion**
   - Extend argparse subcommands; add new options for file formats, chunk size, caching, offline.
   - Implement progress bars (tqdm) and error messaging.
   - Write tests covering CLI parsing and execution with fixtures.
5. **Caching & Offline Mode**
   - Implement local cache for marker DB (using app dirs) with version validation.
   - Allow offline heuristic annotator to run automatically when API key absent.
   - Document environment variables controlling behavior.
6. **Telemetry**
   - Add logging configuration section with examples.
   - Optionally expose `--debug` flag to increase verbosity.
7. **Docs & Notebooks**
   - Write new tutorials (PBMC, integration with pipelines, benchmarking).
   - Update Sphinx/mkdocs (if used) or Markdown docs.
   - Include CLI cheat sheet and troubleshooting FAQ.
8. **Testing & QA**
   - Enhance pytest coverage (use `hypothesis` for edge cases if helpful).
   - Ensure tests run quickly (<5 min) by mocking network or using offline heuristic.
   - Add CLI integration test to GitHub Actions matrix.
9. **Release Automation**
   - Update GitHub Actions to build wheels, run tests, and publish on tag.
   - Sync version bumps across `pyproject.toml`, Docker tags, and docs.
   - Provide post-release validation checklist (install from PyPI, run demo).

## Quality Gates
- pytest, mypy, Ruff, and CLI smoke tests green in CI.
- Wheels install and run on Linux/macOS/Windows (test via `cibuildwheel` or GH Actions).
- Tutorials execute without manual edits (`nbmake` or `pytest --nbmake`).
- Marker DB caching validated in tests with temp directories.
- Offline mode documented and tested.

## Stretch Goals
- Add Airflow/Prefect operators for batch annotation.
- Provide interactive Streamlit dashboards reading AnnData outputs.
- Benchmark suite reporting accuracy vs. GPTCelltype across datasets; publish results.
