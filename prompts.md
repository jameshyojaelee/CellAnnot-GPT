
## P6 – Benchmark Automation CLI & Reports
```text
Objective: automate evaluation against multiple datasets.

Actions:
1. Add `scripts/run_benchmarks.py` reading configs in `evaluation/datasets/*.json`, calling `load_and_run`, saving structured JSON + markdown via `evaluation/report_templates.py`.
2. Store outputs in `docs/reports/YYYYMMDD/benchmark_name.{json,md}`.
3. CI: extend workflow to run benchmarks on schedule (cron) and upload artifacts as build outputs.
4. Include instructions in `docs/benchmarking.md`.
```

## P7 – Deployment Hardening & Secret Templates
```text
Aim: polish release process and secret handling.

Changes:
1. Create `.env.template` listing required vars (OPENAI_API_KEY, CELLANNOT_API_URL, REDIS_URL, ENVIRONMENT).
2. Update `.gitignore` to ensure `.env` stays local.
3. Expand `.github/workflows/ci.yml` to authenticate with GHCR, push image tags, and run Trivy scan.
4. Add `docs/deployment.md` section covering secret sourcing (GH secrets, AWS/GCP secret manager).
```

## P8 – UI Theme & Professional Polish
```text
Goal: make Streamlit app feel like a production analytics tool.

Tune:
1. Configure custom theme (`.streamlit/config.toml`) with consistent palette/typography.
2. Add app-wide navigation header, breadcrumb-style progress indicator, and tooltip icons explaining concepts (validation, confidence).
3. Provide download buttons for results (CSV/JSON/PDF) and screenshot/Export to slides guidance.
4. Incorporate subtle animations (Streamlit elements, spinners) and success/error toast notifications.
```

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
