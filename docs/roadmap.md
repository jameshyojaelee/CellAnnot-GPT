# Roadmap

This roadmap highlights upcoming capabilities and areas of active exploration. Timelines are indicative and prioritised by user feedback and technical feasibility.

## 0–3 Months
- **Expanded marker sources**: ingest immune + tissue-specific datasets (Human Cell Atlas, AZUL) using the existing `scripts/build_marker_db.py` pipeline.
- **UI validation overlays**: expose validation warnings directly on UMAP/heatmap views to speed up triage.
- **CLI ergonomics**: add `gca inspect-report` for summarising JSON outputs and `gca cache clear` to reset bundled assets.
- **Benchmark automation**: nightly `poetry run poe benchmarks` with trend visualisations pushed to `docs/reports/latest/`.

## 3–6 Months
- **Active-learning loops**: surface “mark for review” clusters and feed curator feedback into prompt tuning + knowledge updates.
- **Local LLM backends**: integrate lightweight open models (e.g., `gpt-4o-mini`, `Llama 3.1`) behind the same adapter interface used by `backend/llm/annotator.py`.
- **Ontology alignment tools**: helper CLI to reconcile custom vocabularies with CL IDs, exporting diff reports for curation teams.
- **Pipeline integrations**: official Nextflow and Snakemake modules invoking `gca` commands for reproducible batch jobs.

## 6+ Months
- **Multi-modal support**: extend validation to include protein/ATAC markers when available.
- **Continuous knowledge refresh**: scheduled ingestion of literature-derived markers with citation tracking.
- **Team collaboration features**: shared review queues, annotation comments, and audit logs for regulated environments.

## Recently Shipped
- Offline-aware packaging with bundled assets and the `gca` CLI (see [`docs/install.md`](docs/install.md)).
- Scanpy workflow helpers and CLI parity (`gca scanpy`, [`docs/scanpy_integration.md`](docs/scanpy_integration.md)).
- Enhanced validation guardrails documented in [`docs/operations.md`](docs/operations.md#validation-guardrails).

Suggestions welcome—open an issue or start a discussion with the desired feature, target use case, and timeline.
