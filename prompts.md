
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
