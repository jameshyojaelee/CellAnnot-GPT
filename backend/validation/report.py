"""Utilities to aggregate annotation and validation results into reports."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from typing import Any

from pydantic import BaseModel, Field

from backend.validation.crosscheck import CrosscheckResult
from config.settings import get_settings

FLAG_REASON_MESSAGES: dict[str, str] = {
    "low_marker_overlap": "Low marker overlap with knowledge base markers.",
    "contradictory_markers": "Markers contradict the proposed label.",
    "ontology_mismatch": "Ontology identifier does not match knowledge base records.",
    "label_not_in_kb": "Proposed label not found in marker knowledge base.",
    "species_mismatch": "No marker records for the requested species.",
    "tissue_mismatch": "No marker records for the requested tissue.",
    "no_markers_supplied": "No markers were supplied for validation.",
    "missing_ontology_id": "Ontology identifier missing from annotation payload.",
    "validation_missing": "Validation result missing for this cluster.",
}

CONFIDENCE_ORDER = {"Low": 0, "Medium": 1, "High": 2}


def _calibrate_confidence(original: str | None, support_count: int) -> str:
    settings = get_settings()
    medium_threshold = max(0, settings.confidence_overlap_medium)
    high_threshold = max(medium_threshold, settings.confidence_overlap_high)

    if support_count >= high_threshold > 0:
        computed = "High"
    elif support_count >= medium_threshold > 0:
        computed = "Medium"
    elif support_count > 0:
        computed = "Low"
    else:
        computed = "Low"

    if not original or original.lower() == "unknown":
        return computed

    original_norm = original.title()
    if original_norm not in CONFIDENCE_ORDER:
        return computed

    computed_rank = CONFIDENCE_ORDER[computed]
    original_rank = CONFIDENCE_ORDER[original_norm]
    downgraded_rank = min(original_rank, computed_rank)
    for level, rank in CONFIDENCE_ORDER.items():
        if rank == downgraded_rank:
            return level
    return computed


class ClusterReport(BaseModel):
    cluster_id: str
    annotation: dict[str, Any]
    validation: dict[str, Any] | None = None
    warnings: list[str] = Field(default_factory=list)
    status: str
    confidence: str | None = None


class DatasetSummary(BaseModel):
    total_clusters: int
    supported_clusters: int
    flagged_clusters: int
    unknown_clusters: list[str] = Field(default_factory=list)


class DatasetMetrics(BaseModel):
    support_rate: float
    flagged_rate: float
    unknown_rate: float
    flagged_reasons: dict[str, int] = Field(default_factory=dict)
    confidence_counts: dict[str, int] = Field(default_factory=dict)


class DatasetReport(BaseModel):
    summary: DatasetSummary
    metrics: DatasetMetrics
    clusters: list[ClusterReport]


def _compute_rates(count: int, total: int) -> float:
    return round(count / total, 4) if total else 0.0


def build_structured_report(
    annotations: Iterable[dict[str, Any]],
    crosscheck_results: dict[str, CrosscheckResult],
) -> DatasetReport:
    """Combine annotations with validation findings into a canonical structure."""

    clusters: list[ClusterReport] = []
    supported = 0
    flagged = 0
    unknown_clusters: list[str] = []
    reason_counts: Counter[str] = Counter()
    confidence_counts: Counter[str] = Counter()
    settings = get_settings()

    for annotation in annotations:
        cluster_id = str(annotation.get("cluster_id", "unknown"))
        validation = crosscheck_results.get(cluster_id)
        annotation_data = dict(annotation)
        warnings: list[str] = []

        support_count = validation.support_count if validation else 0
        calibrated_confidence = _calibrate_confidence(
            annotation_data.get("confidence"),
            support_count,
        )

        status = "supported"
        if validation and validation.is_supported:
            supported += 1
        else:
            status = "flagged"
            flagged += 1
            if settings.validation_force_unknown_on_fail:
                if validation and "proposed_label" not in annotation_data:
                    annotation_data["proposed_label"] = validation.primary_label
                annotation_data["primary_label"] = "Unknown or Novel"
                annotation_data["ontology_id"] = None
            calibrated_confidence = "Low"

        annotation_data["confidence"] = calibrated_confidence
        if validation:
            for reason in validation.flag_reasons:
                reason_counts[reason] += 1
                base_message = FLAG_REASON_MESSAGES.get(reason, reason.replace("_", " ").title())
                if base_message:
                    warnings.append(base_message)
                if reason == "contradictory_markers" and validation.contradictory_markers:
                    markers = ", ".join(validation.contradictory_markers.keys())
                    warnings.append(f"Contradictory markers: {markers}")
                if reason == "low_marker_overlap":
                    warnings.append(
                        f"Only {validation.support_count} marker(s) overlap with knowledge base "
                        f"(minimum {max(settings.validation_min_marker_overlap, 1)} required)."
                    )

            if validation.notes:
                for note in validation.notes:
                    if note not in warnings:
                        warnings.append(note)
        else:
            reason_counts["validation_missing"] += 1
            warnings.append(FLAG_REASON_MESSAGES["validation_missing"])
            if settings.validation_force_unknown_on_fail:
                annotation_data["primary_label"] = "Unknown or Novel"
                annotation_data["ontology_id"] = None

        mapping_notes = annotation_data.get("mapping_notes") or annotation_data.get("metadata", {}).get("mapping_notes")
        if mapping_notes:
            annotation_data.setdefault("metadata", {})["mapping_notes"] = mapping_notes
            unresolved = [note for note in mapping_notes if not note.get("target")]
            if unresolved:
                reason_counts["ortholog_unmapped"] += len(unresolved)
                missing = ", ".join(note.get("source", "?") for note in unresolved)
                warnings.append(f"Ortholog mapping missing for: {missing}")

        # Deduplicate warnings while preserving order.
        deduped_warnings: list[str] = []
        seen_warnings: set[str] = set()
        for warning in warnings:
            if warning and warning not in seen_warnings:
                deduped_warnings.append(warning)
                seen_warnings.add(warning)

        if annotation_data.get("primary_label") in (None, "", "Unknown", "Unknown or Novel"):
            unknown_clusters.append(cluster_id)

        confidence_label = (annotation_data.get("confidence") or "Unknown").title()
        annotation_data["confidence"] = confidence_label
        confidence_counts[confidence_label] += 1

        clusters.append(
            ClusterReport(
                cluster_id=cluster_id,
                annotation=annotation_data,
                validation=validation.to_dict() if validation else None,
                warnings=deduped_warnings,
                status=status,
                confidence=annotation_data["confidence"],
            )
        )

    total = len(clusters)
    summary = DatasetSummary(
        total_clusters=total,
        supported_clusters=supported,
        flagged_clusters=flagged,
        unknown_clusters=unknown_clusters,
    )
    metrics = DatasetMetrics(
        support_rate=_compute_rates(supported, total),
        flagged_rate=_compute_rates(flagged, total),
        unknown_rate=_compute_rates(len(unknown_clusters), total),
        flagged_reasons=dict(reason_counts),
        confidence_counts=dict(confidence_counts),
    )
    return DatasetReport(summary=summary, metrics=metrics, clusters=clusters)


def render_text_report(report: DatasetReport | dict[str, Any]) -> str:
    """Render a human-readable multi-line summary from the structured report."""

    dataset = report if isinstance(report, DatasetReport) else DatasetReport.model_validate(report)

    summary = dataset.summary
    lines = [
        "GPT Cell Annotator Validation Report",
        "================================",
        f"Total clusters: {summary.total_clusters}",
        f"Supported: {summary.supported_clusters}",
        f"Flagged: {summary.flagged_clusters}",
    ]

    if summary.unknown_clusters:
        lines.append(f"Unknown clusters: {', '.join(summary.unknown_clusters)}")

    for cluster in dataset.clusters:
        if not cluster.warnings:
            continue
        lines.append(f"- Cluster {cluster.cluster_id} ({cluster.status}):")
        for warn in cluster.warnings:
            lines.append(f"    - {warn}")

    if len(lines) == 4 and not summary.unknown_clusters:
        lines.append("All clusters supported by marker knowledge base.")

    return "\n".join(lines)


__all__ = [
    "ClusterReport",
    "DatasetSummary",
    "DatasetMetrics",
    "DatasetReport",
    "build_structured_report",
    "render_text_report",
]
