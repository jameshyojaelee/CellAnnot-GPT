# Operations & Observability

## Logging
- CellAnnot-GPT uses `structlog` to emit JSON logs enriched with `trace_id`, `method`, `path`, `status_code`, `duration_ms`. Each response includes `X-Trace-Id` and `X-Process-Time-ms` headers for correlation.
- Tuning: set `LOG_LEVEL` (INFO, DEBUG, etc.) in environment configuration.

## Log Aggregation
- **Local**: pipe logs to a file or use `docker compose logs` for rapid debugging.
- **ELK / OpenSearch**: ship container stdout to Logstash/Fluentbit. Index by `trace_id` to reconstruct request/response flows.
- **Cloud providers**: leverage AWS CloudWatch, GCP Cloud Logging, or Azure Monitor. Map JSON logs into structured dashboards and alerts.

## Tracing & Metrics
- `trace_id` header enables propagation across upstream services. Extend by forwarding to downstream APIs and central traces (e.g., OpenTelemetry).
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

