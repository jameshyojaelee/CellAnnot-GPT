"""Utilities to aggregate annotation and validation results into reports."""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, Iterable, List, Optional, Union

from pydantic import BaseModel, Field

from backend.validation.crosscheck import CrosscheckResult


class ClusterReport(BaseModel):
    cluster_id: str
    annotation: Dict[str, Any]
    validation: Optional[Dict[str, Any]] = None
    warnings: List[str] = Field(default_factory=list)
    status: str
    confidence: Optional[str] = None


class DatasetSummary(BaseModel):
    total_clusters: int
    supported_clusters: int
    flagged_clusters: int
    unknown_clusters: List[str] = Field(default_factory=list)


class DatasetMetrics(BaseModel):
    support_rate: float
    flagged_rate: float
    unknown_rate: float
    flagged_reasons: Dict[str, int] = Field(default_factory=dict)
    confidence_counts: Dict[str, int] = Field(default_factory=dict)


class DatasetReport(BaseModel):
    summary: DatasetSummary
    metrics: DatasetMetrics
    clusters: List[ClusterReport]


def _compute_rates(count: int, total: int) -> float:
    return round(count / total, 4) if total else 0.0


def build_structured_report(
    annotations: Iterable[Dict[str, Any]],
    crosscheck_results: Dict[str, CrosscheckResult],
) -> DatasetReport:
    """Combine annotations with validation findings into a canonical structure."""

    clusters: List[ClusterReport] = []
    supported = 0
    flagged = 0
    unknown_clusters: List[str] = []
    reason_counts: Counter[str] = Counter()
    confidence_counts: Counter[str] = Counter()

    for annotation in annotations:
        cluster_id = str(annotation.get("cluster_id", "unknown"))
        validation = crosscheck_results.get(cluster_id)
        warnings: List[str] = []
        confidence = (annotation.get("confidence") or "Unknown").title()
        confidence_counts[confidence] += 1

        status = "supported"
        if validation:
            if validation.is_supported:
                supported += 1
            else:
                status = "flagged"
                flagged += 1
            if validation.ontology_mismatch:
                warnings.append("Ontology mismatch between annotation and DB")
                reason_counts["ontology_mismatch"] += 1
            if validation.contradictory_markers:
                markers = ", ".join(validation.contradictory_markers.keys())
                warnings.append(f"Contradictory markers: {markers}")
                reason_counts["contradictory_markers"] += 1
            if validation.notes:
                warnings.extend(validation.notes)
                reason_counts["notes"] += len(validation.notes)
        else:
            status = "flagged"
            flagged += 1
            warnings.append("Validation result missing")
            reason_counts["validation_missing"] += 1

        if annotation.get("primary_label") in (None, "", "Unknown", "Unknown or Novel"):
            unknown_clusters.append(cluster_id)

        clusters.append(
            ClusterReport(
                cluster_id=cluster_id,
                annotation=annotation,
                validation=validation.to_dict() if validation else None,
                warnings=warnings,
                status=status,
                confidence=confidence,
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


def render_text_report(report: Union[DatasetReport, Dict[str, Any]]) -> str:
    """Render a human-readable multi-line summary from the structured report."""

    dataset = report if isinstance(report, DatasetReport) else DatasetReport.model_validate(report)

    summary = dataset.summary
    lines = [
        "CellAnnot-GPT Validation Report",
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
