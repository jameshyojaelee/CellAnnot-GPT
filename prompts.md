

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
