# Mega Prompt: Build GPT Cell Annotator R/Seurat Companion Package

You are GPT-5 Codex acting as a senior engineer tasked with shipping a production-ready R package that connects Seurat workflows to GPT Cell Annotator. Deliver the package end-to-end with robust tests, docs, and release automation.

## Goals
- Expose a user-friendly R API (`gptca_annotate_*`) that accepts Seurat `FindAllMarkers()` tables or custom gene lists and returns validated labels plus confidence/status metadata.
- Keep all core business logic inside the existing Python backend; the R package should communicate via HTTPS with the FastAPI service (preferred) but support an offline fallback by shelling out to the CLI if needed.
- Surface validation guardrails (downgrades, warnings) in Seurat metadata and provide plotting helpers to visualize annotations on UMAP.
- Ship complete documentation (pkgdown site + vignettes), CI, and release automation.

## Constraints / Assumptions
- Repository: `github.com/.../GPT-cell-annotator` (Python). The new R package lives in `clients/r`.
- Backend REST endpoints: `/annotate_cluster` and `/annotate_batch`, request/response schema defined in `schemas/annotation.py`.
- Prefer HTTPS calls using `httr2`. Do not reimplement the LLM logic in R.
- Package must support R >= 4.1. Provide fallbacks or meaningful errors for Windows/macOS (Rtools/clang guidance).
- Provide knobs for `species`, `tissue`, marker limits, timeout, retry, and offline CLI fallback (`gca annotate`).
- Tests must be runnable on CI without real OpenAI keys: mock REST responses and record fixtures.
- Document costs and rate-limit considerations; warn users not to leak API keys.

## Deliverables
1. R package structure in `clients/r/gptcellannotator` with:
   - `DESCRIPTION`, `NAMESPACE`, `R/` sources, `tests/`, `man/`, `vignettes/`.
   - Core functions: `gptca_annotate_markers`, `gptca_annotate_seurat`, `gptca_add_metadata`, `gptca_plot_umap`.
   - Internal helpers for schema validation, REST/CLI transport, and configuration (`GptcaConfig` S3 class).
2. Automated test suite:
   - Unit tests mocking HTTP responses (`httptest2`).
   - Integration test using bundled PBMC markers and mocked backend.
   - Lint/style checks (`lintr`, `styler`).
3. Documentation:
   - Roxygen comments for all exports.
   - Vignette: “Annotate Seurat Clusters with GPT Cell Annotator”.
   - pkgdown configuration and GitHub Actions workflow.
4. CLI/offline bridge leveraging Python CLI: detect CLI availability, run `gca annotate` on temp CSV, parse JSON, map fields back to R objects.
5. Developer ergonomics:
   - Makefile or `devtools::` wrappers.
   - CI workflow(s) under `.github/workflows` covering R CMD check, pkgdown deploy, unit tests.
6. Release guide in `clients/r/RELEASING.md` describing versioning, dependency pinning, and publishing to CRAN (if feasible) or GitHub releases.

## Execution Steps
1. **Discovery**
   - Inspect Python schema definitions (`schemas/annotation.py`, `docs/api_reference.md`) and confirm field naming.
   - Gather existing demo markers (`data/demo/pbmc_markers.csv`) and expected outputs from tests.
2. **Scaffold Package**
   - Use `usethis::create_package("clients/r/gptcellannotator")`.
   - Configure package options, dependencies (`httr2`, `jsonlite`, `seurat`, `cli`, `rlang`, `purrr`, `ggplot2`).
   - Set up `DESCRIPTION`, licensing, `NAMESPACE` exports.
3. **Implement Configuration Layer**
   - Provide `gptca_config()` to store base URL, API key, CLI path, retries, timeouts.
   - Support environment variable overrides (`GPTCA_BASE_URL`, `GPTCA_API_KEY`).
4. **Transport + Schema**
   - Implement HTTP client with `httr2` (JSON encode/decode, error handling, exponential backoff).
   - Define response parsing that maps to tidy R structures (`tibble`).
   - CLI fallback: write markers to temp CSV, call `system2("gca", args...)`, parse output.
5. **Seurat Integration**
   - Convert Seurat marker tables (`FindAllMarkers`) to required payload (top-N genes per cluster).
   - Add helper to merge results into `obj@meta.data`, preserving validation warnings and ontology IDs.
   - Provide `gptca_plot_umap(obj, label_col="gptca_label")` wrapper using `DimPlot`.
6. **Validation & Edge Cases**
   - Handle missing markers, empty clusters, large datasets (chunk requests).
   - Make guardrails user-visible: warnings column, unknown downgrades, rationale text.
   - Include functions to inspect raw JSON report.
7. **Testing**
   - Unit tests for payload construction, HTTP/CLI transports, parsing.
   - Snapshot tests for typical responses.
   - Integration test running through Seurat PBMC example with mocked backend.
8. **Docs & Examples**
   - Roxygen documentation.
   - Vignette walking through full pipeline (setup, annotate, visualize, offline fallback).
   - Troubleshooting section (rate limits, HTTP errors).
9. **Automation**
   - GitHub Actions workflow for R CMD check (`macOS-latest`, `ubuntu-latest`).
   - pkgdown site build and GitHub Pages deployment.
   - Pre-commit hooks (optional) for formatting.
10. **Release Preparation**
   - Version everything, update changelog, create release checklist.
   - Coordinate with Python team on API compatibility guarantees.

## Quality Gates
- All R CMD check steps pass with `--as-cran`.
- Test coverage >= 80% for exported functions (measure via `covr`).
- Vignette renders without manual intervention.
- CLI fallback tested on CI (mock CLI binary if needed).
- Documented manual QA steps (Seurat PBMC run against staging backend).

## Stretch Goals
- Optional: wrap backend benchmark endpoints and produce comparative plots.
- Optional: expose shiny gadget for interactive annotation review.
- Optional: publish to CRAN after internal review.
