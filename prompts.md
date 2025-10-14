# GPT-5-Codex Prompt Playbook · CellAnnot-GPT

Each section is a ready-to-run instruction block for Codex. Run them one at a time, reviewing diffs and tests before proceeding.



## 2. LLM Orchestration & Schema Hardening

**Prompt**

You are refactoring the CellAnnot-GPT annotator.
Objectives:
1. Remove the duplicate system message in `backend/llm/annotator.py` and introduce structured logging around request/response payloads (without leaking secrets).
2. Switch OpenAI calls to use `response_format={"type":"json_object"}` (or function calling) tied to `schemas/annotation_result.schema.json`.
3. Implement a JSON Schema validator so malformed outputs trigger the `MockAnnotator` fallback with an appended warning.
4. Extend prompts in `backend/llm/prompts.py` so the model must cite at least one knowledge-base evidence string.

Remember to update tests in `tests/test_annotator.py` and add new cases for schema validation. Run `poetry run pytest tests/test_annotator.py` after coding.

---

## 3. API & Cache Reliability

**Prompt**

Focus on the FastAPI layer.
Tasks:
1. In `backend/api/main.py`, honor `return_validated=False` by returning the raw annotation payload (and reserve crosscheck only for the validated path).
2. Ensure cache keys encode whether validation was requested; avoid mixing validated and non-validated reports.
3. Add structured error handling so Redis outages degrade to a warning log and continue without caching.
4. Extend `backend/cache/async_cache.py` to accept async callables directly.
5. Update `tests/test_api.py` and `tests/test_cache.py` to cover the new semantics.

Run `poetry run pytest tests/test_api.py tests/test_cache.py` before finishing.

---

## 4. Benchmark Guardrails

**Prompt**

Make benchmarks a first-class signal.
Deliverables:
1. Accept multiple datasets via globbing in `scripts/run_benchmarks.py` and output accuracy deltas versus the previous snapshot stored in `docs/reports/latest`.
2. Emit markdown that highlights regressions >5% in bold.
3. Update `.github/workflows/ci.yml` so scheduled runs upload both raw JSON and diff summaries, and fail the job on regression unless `ALLOW_BENCHMARK_REGRESSIONS` is true.
4. Create a helper in `evaluation/report_templates.py` that renders sparkline-ready CSV for dashboards.

Add or update tests under `tests/test_benchmark_runner.py`. Run `poetry run pytest tests/test_benchmark_runner.py`.

---

## 5. Signature Streamlit Experience

**Prompt**

Elevate the UI.
Goals:
1. In `frontend/streamlit_app.py`, add a comparison mode to overlay two batches and surface validation contradictions visually (Altair or Plotly is fine).
2. Highlight markers that triggered warnings, linking back to the knowledge base source (from the API payload).
3. Provide call-to-action panels (“Next experiments”, “Markers to validate”) powered by validation notes.
4. Extract shared formatting helpers into `frontend/utils.py` with unit tests in `tests/test_frontend.py`.
5. Update README screenshots or add a short Loom link in `docs/demo.md`.

End with `poetry run pytest tests/test_frontend.py` and, if feasible, `poetry run streamlit run frontend/streamlit_app.py --server.headless true --browser.gatherUsageStats false`.

---
