"""Utilities for Streamlit frontend to interact with backend/API."""

from __future__ import annotations

import functools
import os
from collections.abc import Sequence
from typing import Any

import requests

API_BASE_URL = os.environ.get("CELLANNOT_API_URL", "http://127.0.0.1:8000")


@functools.lru_cache(maxsize=32)
def get_health() -> dict[str, Any]:
    """Check API health endpoint."""

    response = requests.get(f"{API_BASE_URL}/health", timeout=10)
    response.raise_for_status()
    return response.json()


def annotate_cluster(payload: dict[str, Any]) -> dict[str, Any]:
    """Call backend for single cluster annotation."""

    response = requests.post(
        f"{API_BASE_URL}/annotate_cluster",
        json=payload,
        timeout=60,
    )
    response.raise_for_status()
    return response.json()


def annotate_batch(payload: dict[str, Any]) -> dict[str, Any]:
    """Call backend for batch annotation."""

    response = requests.post(
        f"{API_BASE_URL}/annotate_batch",
        json=payload,
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


def annotate_cluster_api(payload: dict[str, Any]) -> dict[str, Any]:
    response = requests.post(
        f"{API_BASE_URL}/annotate_cluster",
        json=payload,
        timeout=60,
    )
    response.raise_for_status()
    return response.json()


def status_badge(label: str, status: str) -> str:
    colors = {
        "ok": "#3CB371",
        "warn": "#FFA500",
        "error": "#DC143C",
    }
    color = colors.get(status, "#808080")
    return (
        f"<span style='padding:4px 8px;border-radius:12px;background:{color};color:white;'>"
        f"{label}</span>"
    )


def format_summary(report: dict[str, Any]) -> str:
    summary = report.get("summary", {})
    metrics = report.get("metrics", {})
    support_rate = metrics.get("support_rate")
    flagged_rate = metrics.get("flagged_rate")
    unknown_rate = metrics.get("unknown_rate")
    total = summary.get("total_clusters", 0)
    supported = summary.get("supported_clusters", 0)
    flagged = summary.get("flagged_clusters", 0)
    unknown = len(summary.get("unknown_clusters", []))

    def pct(value: float | None) -> str:
        return f"{value * 100:.1f}%" if value is not None else "n/a"

    return (
        f"Total: {total} | "
        f"Supported: {supported} ({pct(support_rate)}) | "
        f"Flagged: {flagged} ({pct(flagged_rate)}) | "
        f"Unknown: {unknown} ({pct(unknown_rate)})"
    )


def knowledge_base_url(marker: str, source: str | None = None) -> str:
    """Heuristic mapping from marker to a public knowledge base URL."""

    symbol = marker.strip().upper()
    if source and "cellmarker" in source.lower():
        return f"https://biocc.hrbmu.edu.cn/CellMarker/search.jsp?g={symbol}"
    if source and "panglao" in source.lower():
        return f"https://panglaodb.se/markers/{symbol}"
    return f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={symbol}"


def format_marker_links(markers: Sequence[str], source: str | None = None) -> str:
    """Return Markdown string with marker hyperlinks."""

    unique: list[str] = []
    seen: set[str] = set()
    for marker in markers:
        marker = marker.strip()
        if not marker:
            continue
        if marker.upper() in seen:
            continue
        seen.add(marker.upper())
        unique.append(marker)
    if not unique:
        return "—"
    links = [f"[{m}]({knowledge_base_url(m, source)})" for m in unique]
    return ", ".join(links)


def collect_marker_highlights(cluster: dict[str, Any]) -> dict[str, list[str]]:
    """Extract marker lists used for warnings and support visuals."""

    validation = cluster.get("validation") or {}
    contradictory = validation.get("contradictory_markers") or {}
    missing = validation.get("missing_markers") or []
    supporting = validation.get("supporting_markers") or []

    warning_markers: list[str] = sorted({*contradictory.keys(), *missing})
    supporting_markers: list[str] = sorted(set(supporting))
    return {
        "warning_markers": warning_markers,
        "supporting_markers": supporting_markers,
        "contradictory_map": contradictory,
        "missing_markers": missing,
    }


def build_call_to_action(cluster: dict[str, Any]) -> dict[str, list[str]]:
    """Generate follow-up suggestions based on validation output."""

    validation = cluster.get("validation") or {}
    warnings = cluster.get("warnings") or []
    contradictory = validation.get("contradictory_markers") or {}
    missing = validation.get("missing_markers") or []
    confidence = (cluster.get("confidence") or "").lower()

    experiments: list[str] = []
    markers_to_validate: list[str] = []

    if contradictory:
        markers_to_validate.extend(contradictory.keys())
        experiments.append(
            "Design flow cytometry or staining panel to resolve conflicting markers."
        )
    if missing:
        markers_to_validate.extend(missing)
        experiments.append(
            "Collect additional markers or deeper sequencing for low-evidence genes."
        )
    if warnings and "novel" in " ".join(warnings).lower():
        experiments.append(
            "Review novel cluster with domain experts before assigning a definitive label."
        )
    if confidence in {"low", "unknown"} and not experiments:
        experiments.append(
            "Plan orthogonal assay (e.g., CITE-seq) to firm up low-confidence calls."
        )

    if not experiments:
        experiments.append("Document validated markers and archive the run for traceability.")

    dedup_markers = sorted({m.upper(): m for m in markers_to_validate}.values())
    return {
        "next_experiments": sorted(set(experiments)),
        "markers_to_validate": dedup_markers,
    }


__all__ = [
    "annotate_batch",
    "annotate_cluster",
    "annotate_cluster_api",
    "build_call_to_action",
    "collect_marker_highlights",
    "format_marker_links",
    "format_summary",
    "get_health",
    "knowledge_base_url",
    "status_badge",
]
