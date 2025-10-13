"""Utilities to aggregate annotation and validation results into reports."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List

from backend.validation.crosscheck import CrosscheckResult


def build_structured_report(
    annotations: Iterable[Dict[str, Any]],
    crosscheck_results: Dict[str, CrosscheckResult],
) -> Dict[str, Any]:
    """Combine annotations with validation findings into a canonical structure."""

    clusters_payload = []
    supported = 0
    flagged = 0
    unknown_clusters: List[str] = []

    for annotation in annotations:
        cluster_id = str(annotation.get("cluster_id", "unknown"))
        validation = crosscheck_results.get(cluster_id)
        warnings: List[str] = []

        if validation:
            if validation.is_supported:
                supported += 1
            else:
                flagged += 1
            if validation.ontology_mismatch:
                warnings.append("Ontology mismatch between annotation and DB")
            if validation.contradictory_markers:
                markers = ", ".join(validation.contradictory_markers.keys())
                warnings.append(f"Contradictory markers: {markers}")
            if validation.notes:
                warnings.extend(validation.notes)
        else:
            flagged += 1
            warnings.append("Validation result missing")

        if annotation.get("primary_label") in (None, "", "Unknown", "Unknown or Novel"):
            unknown_clusters.append(cluster_id)

        clusters_payload.append(
            {
                "cluster_id": cluster_id,
                "annotation": annotation,
                "validation": validation.to_dict() if validation else None,
                "warnings": warnings,
            }
        )

    summary = {
        "total_clusters": len(clusters_payload),
        "supported_clusters": supported,
        "flagged_clusters": flagged,
        "unknown_clusters": unknown_clusters,
    }

    return {
        "summary": summary,
        "clusters": clusters_payload,
    }


def render_text_report(structured_report: Dict[str, Any]) -> str:
    """Render a human-readable multi-line summary from the structured report."""

    summary = structured_report.get("summary", {})
    clusters: List[Dict[str, Any]] = structured_report.get("clusters", [])

    lines = [
        "CellAnnot-GPT Validation Report",
        "================================",
        f"Total clusters: {summary.get('total_clusters', 0)}",
        f"Supported: {summary.get('supported_clusters', 0)}",
        f"Flagged: {summary.get('flagged_clusters', 0)}",
    ]

    unknown = summary.get("unknown_clusters") or []
    if unknown:
        lines.append(f"Unknown clusters: {', '.join(unknown)}")

    for cluster in clusters:
        warnings = cluster.get("warnings") or []
        if not warnings:
            continue
        cluster_id = cluster.get("cluster_id")
        lines.append(f"- Cluster {cluster_id}:")
        for warn in warnings:
            lines.append(f"    - {warn}")

    if len(lines) == 4 and not unknown:
        lines.append("All clusters supported by marker knowledge base.")

    return "\n".join(lines)


__all__ = ["build_structured_report", "render_text_report"]
