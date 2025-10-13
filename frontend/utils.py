"""Utilities for Streamlit frontend to interact with backend/API."""

from __future__ import annotations

import functools
import os
from typing import Any, Dict, List

import requests
import streamlit as st


def _get_api_base_url() -> str:
    try:
        return st.secrets["api_base_url"]  # type: ignore[index]
    except Exception:  # noqa: BLE001 - secrets may be absent in tests
        return os.environ.get("CELLANNOT_API_URL", "http://127.0.0.1:8000")


API_BASE_URL = _get_api_base_url()


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


def status_badge(label: str, status: str) -> str:
    colors = {
        "ok": "#3CB371",
        "warn": "#FFA500",
        "error": "#DC143C",
    }
    color = colors.get(status, "#808080")
    return f"<span style='padding:4px 8px;border-radius:12px;background:{color};color:white;'>{label}</span>"


def format_summary(result: Dict[str, Any]) -> str:
    summary = result.get("summary", {})
    return (
        f"Total: {summary.get('total_clusters', 0)} | "
        f"Supported: {summary.get('supported_clusters', 0)} | "
        f"Flagged: {summary.get('flagged_clusters', 0)} | "
        f"Unknown: {len(summary.get('unknown_clusters', []))}"
    )
