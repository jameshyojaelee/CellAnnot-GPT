"""Utilities for Streamlit frontend to interact with backend/API."""

from __future__ import annotations

import functools
import os
from typing import Any, Dict, List

import requests


API_BASE_URL = os.environ.get("CELLANNOT_API_URL", "http://127.0.0.1:8000")


@functools.lru_cache(maxsize=32)
def get_health() -> Dict[str, Any]:
    """Check API health endpoint."""

    response = requests.get(f"{API_BASE_URL}/health", timeout=10)
    response.raise_for_status()
    return response.json()


def annotate_cluster(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Call backend for single cluster annotation."""

    response = requests.post(
        f"{API_BASE_URL}/annotate_cluster",
        json=payload,
        timeout=60,
    )
    response.raise_for_status()
    return response.json()


def annotate_batch(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Call backend for batch annotation."""

    response = requests.post(
        f"{API_BASE_URL}/annotate_batch",
        json=payload,
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


def annotate_cluster_api(payload: Dict[str, Any]) -> Dict[str, Any]:
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
    return f"<span style='padding:4px 8px;border-radius:12px;background:{color};color:white;'>{label}</span>"


def format_summary(report: Dict[str, Any]) -> str:
    summary = report.get("summary", {})
    metrics = report.get("metrics", {})
    support_rate = metrics.get("support_rate")
    flagged_rate = metrics.get("flagged_rate")
    unknown_rate = metrics.get("unknown_rate")
    total = summary.get("total_clusters", 0)
    supported = summary.get("supported_clusters", 0)
    flagged = summary.get("flagged_clusters", 0)
    unknown = len(summary.get("unknown_clusters", []))

    def pct(value: Optional[float]) -> str:
        return f"{value * 100:.1f}%" if value is not None else "n/a"

    return (
        f"Total: {total} | "
        f"Supported: {supported} ({pct(support_rate)}) | "
        f"Flagged: {flagged} ({pct(flagged_rate)}) | "
        f"Unknown: {unknown} ({pct(unknown_rate)})"
    )
