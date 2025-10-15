# Operations & Observability

## Logging
- GPT Cell Annotator uses `structlog` to emit JSON logs enriched with `trace_id`, `method`, `path`, `status_code`, `duration_ms`. Each response includes `X-Trace-Id` and `X-Process-Time-ms` headers for correlation.
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

## Quality & Validation
- **Marker overlap threshold:** set `VALIDATION_MIN_MARKER_OVERLAP` (defaults to 2) to require a minimum number of shared markers before a label is considered supported. Failing clusters are automatically downgraded to `Unknown or Novel` when `VALIDATION_FORCE_UNKNOWN_ON_FAIL` is left enabled.
- **Confidence calibration:** confidence bands are recalculated from marker overlap counts using `CONFIDENCE_OVERLAP_MEDIUM` / `CONFIDENCE_OVERLAP_HIGH`. Overly optimistic LLM scores are downgraded based on these thresholds.
- **Flag reasons:** validation output includes machine- and human-readable reasons such as `low_marker_overlap`, `species_mismatch`, and `missing_ontology_id`. The UI surfaces these reasons and metrics under “Flagged reasons”.
- **Original suggestions:** when a cluster is flagged, the report retains `annotation.proposed_label` alongside the downgraded `primary_label` so reviewers can still inspect the model’s initial suggestion.
- **Schema enforcement:** the LLM response schema requires ontology identifiers and marker lists. Validation failures issue structured logs containing trace IDs for easier troubleshooting.
