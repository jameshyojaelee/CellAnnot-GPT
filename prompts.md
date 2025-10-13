
## P9 – Observability & Structured Logging
```text
Goal: increase insight into system behavior.

Implement:
1. Adopt `structlog` or `loguru` in backend to emit JSON logs (request id, llm_mode, cache hits).
2. Add request tracing middleware storing trace ids in responses.
3. Document log aggregation strategy (e.g., shipping to ELK) in `docs/operations.md`.
4. Extend tests validating log output format and failure handling.
```

## P10 – Documentation & Onboarding Refresh
```text
Purpose: help new contributors/users ramp quickly.

Tasks:
1. Create `docs/getting_started.md` covering env setup, running tests, sample workflows.
2. Add architecture diagram (Mermaid or draw.io) in `docs/architecture.md`.
3. Generate API reference via FastAPI docs export (e.g., `docs/api_reference.md`).
4. Link docs from README badges and Streamlit sidebar help modal.
```

---

### Usage Notes
- Before running each prompt, sync repo and ensure tests pass.
- After GPT‑5 responds, apply patches carefully, run `pytest`, `ruff`, and `mypy` as needed.
- Update this file if prompts evolve so the playbook stays current.
