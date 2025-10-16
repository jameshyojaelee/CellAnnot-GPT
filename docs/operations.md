# Operations & Observability

## Logging
- GPT Cell Annotator uses `structlog` to emit JSON logs enriched with `trace_id`, `method`, `path`, `status_code`, `duration_ms`. Each response includes `X-Trace-Id` and `X-Process-Time-ms` headers for correlation.
- Tuning: set `LOG_LEVEL` (INFO, DEBUG, etc.) in environment configuration.

## Offline Assets & Cache
- `GCA_MARKER_DB_PATH` overrides the auto-materialised marker database path (useful when mounting corporate KBs or NFS volumes).
- `GCA_CACHE_DIR` (and the CLI `--cache-dir` flag) enables the disk-backed annotation cache; reuse annotations across repeated reruns or air-gapped environments.
- Cache directories contain JSON blobs keyed by cluster payloads. Clean the directory to invalidate stale responses.

## Log Aggregation
- **Local**: pipe logs to a file or use `docker compose logs` for rapid debugging.
- **ELK / OpenSearch**: ship container stdout to Logstash/Fluentbit. Index by `trace_id` to reconstruct request/response flows.
- **Cloud providers**: leverage AWS CloudWatch, GCP Cloud Logging, or Azure Monitor. Map JSON logs into structured dashboards and alerts.

## Tracing & Metrics
- `trace_id` header enables propagation across upstream services. Extend by forwarding to downstream APIs and central traces (e.g., OpenTelemetry).
- Install the `[api]` extra to expose default Prometheus metrics (`gca_scanpy_batches_total`, `gca_scanpy_clusters_total`, `gca_scanpy_batch_duration_seconds`). Scrape them alongside existing API counters.
- Capture additional metrics (success rate, latency percentiles) via Prometheus exporters or APM tools.

## Alerts & Playbooks
- Set thresholds on error rates, annotation latency, or cache miss ratios. Combine logs + metrics for actionable alerts.
- Maintain runbooks describing remediation steps (restart pods, rotate keys, check Redis health).

## Data Retention
- Configure log retention based on compliance (typically 30–90 days). Ensure personally identifiable data is not logged; marker lists are not sensitive but dataset names may be.

## Observability Toolkit Roadmap
- Integrate OpenTelemetry SDK for distributed traces.
- Add SLO dashboards (Grafana) fed by metrics and logs.
- Automate synthetic probes to verify health endpoints.

## Validation Guardrails

Validation runs after every annotation (CLI, API, notebooks) and governs which labels are promoted, downgraded, or flagged for review.

### Core thresholds

| Setting | Default | Purpose |
| --- | --- | --- |
| `VALIDATION_MIN_MARKER_OVERLAP` | `2` | Minimum shared markers between prediction and marker DB to remain “supported”. |
| `VALIDATION_FORCE_UNKNOWN_ON_FAIL` | `true` | Downgrades clusters with insufficient evidence to `Unknown or Novel`. |
| `CONFIDENCE_OVERLAP_MEDIUM` / `CONFIDENCE_OVERLAP_HIGH` | `2` / `3` | Translate overlap counts into `Low`, `Medium`, `High` confidence bands. |

Set these in `.env` (loaded by `config/settings.py`) or export them before starting the API/CLI (`VALIDATION_MIN_MARKER_OVERLAP=3 gca annotate ...`).

### Status & warnings

- `supported`: Validation thresholds met; status is reflected in `_status` columns and reports.
- `flagged_low_support`, `flagged_species_mismatch`, `flagged_missing_ontology`, etc.: Stored in `annotation.warnings` and surfaced by the Streamlit UI.
- `annotation.proposed_label`: Original LLM suggestion retained even when the primary label is downgraded.

Structured logs contain the `trace_id`, `cluster_id`, and warning codes—ideal for alerting or audit pipelines.

### Retrieval and synonym controls

- `RAG_ENABLED` / `RAG_TOP_K` / `RAG_MIN_OVERLAP` tune the retrieval candidates added to prompts.
- `SYNONYM_CONFIG_PATH` and `ORTHOLOG_MAPPING_PATH` point to custom synonym/ortholog files.
- `SYNONYM_ENABLE_ORTHOLOGS=false` disables cross-species expansion for debugging.

Monitor the effect of these settings via `docs/benchmarks.md` reports and the CLI smoke tests in `tests/test_cli.py`.
