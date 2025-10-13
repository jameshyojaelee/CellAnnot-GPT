# Architecture Overview

## Diagram (Textual)

```
User UI (Streamlit)
      │
      ▼
FastAPI Service ──▶ Annotation Engine (LLM Client)
      │                      │
      │                      ├── Prompt Templates
      │                      └── Rate-limited OpenAI Interface
      │
      ├── Validation Layer ──▶ Marker Knowledge Store (SQLite/Parquet)
      │                              ▲
      │                              │
      └── Reporting Output ◀─────────┘
Benchmark Toolkit ↔ FastAPI/Engine (offline evaluation)
```

## Components

- **Frontend (Streamlit):** Collects marker inputs, invokes API/local engine, renders annotations, explanations, and warnings.
- **API Service (FastAPI):** Exposes REST endpoints for single/batch annotation, health checks, and integrates logging, auth, and validation middleware.
- **Annotation Engine:** Encapsulates prompt construction, OpenAI GPT-4o calls, retry/backoff logic, and response parsing.
- **Validation Layer:** Cross-checks LLM predictions with the marker knowledge store, flags contradictions, and produces structured reports.
- **Marker Knowledge Store:** Normalized SQLite/Parquet database built from PanglaoDB, CellMarker, and curated literature entries.
- **Benchmark Toolkit:** Runs evaluation datasets through the engine, computes metrics, and generates markdown/HTML summaries.
- **Config & Secrets:** Centralized via `config/settings.py` using environment variables for API keys, model selection, caching options.

## Data Flow

1. User uploads marker lists or selects clusters via UI/API.
2. API forwards payload to Annotation Engine, which formats prompts and queries GPT-4o.
3. Engine returns candidate labels with rationales; Validation Layer verifies against marker DB.
4. Responses (including warnings and confidence) return to API/UI and are optionally logged for auditing.
5. Benchmark Toolkit reuses the engine to score labeled datasets and informs prompt/DB updates.

## Non-Functional Considerations

- **Scalability:** Stateless API with optional async queue for batch jobs; caching frequent annotations.
- **Security:** Marker-only inputs, no PII; secrets via env vars; audit logs configurable.
- **Extensibility:** LLM client abstracted to support local models; knowledge store loaders modular for new sources.
