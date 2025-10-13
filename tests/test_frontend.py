"""Lightweight frontend utility tests.

Streamlit does not ship with an official headless testing harness, so we
exercise the helper functions (API adapters, formatting helpers) directly.
If streamlit-testing-tools is added later, these tests can be expanded to
snapshot page layouts.
"""

from __future__ import annotations

import json

import pytest

from frontend import utils


class DummyResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP error {self.status_code}")

    def json(self):
        return self._payload

    @property
    def text(self):
        return json.dumps(self._payload)


def test_format_summary_builds_human_string():
    summary = utils.format_summary(
        {
            "summary": {
                "total_clusters": 3,
                "supported_clusters": 2,
                "flagged_clusters": 1,
                "unknown_clusters": ["1"],
            },
            "metrics": {
                "support_rate": 2 / 3,
                "flagged_rate": 1 / 3,
                "unknown_rate": 1 / 3,
            },
        }
    )
    assert "Total: 3" in summary
    assert "Flagged: 1 (33.3%" in summary


def test_status_badge_returns_html():
    html = utils.status_badge("API", "ok")
    assert "API" in html
    assert "background" in html


def test_api_helpers(monkeypatch):
    calls = {}

    def fake_get(url, timeout=10):
        calls["get"] = (url, timeout)
        return DummyResponse({"status": "ok", "llm_mode": "mock", "cache_enabled": False})

    def fake_post(url, data=None, headers=None, timeout=None, json=None, **kwargs):
        payload = json if json is not None else data
        calls.setdefault("post", []).append((url, payload, headers, timeout))
        return DummyResponse({"result": {}})

    monkeypatch.setattr(utils.requests, "get", fake_get)
    monkeypatch.setattr(utils.requests, "post", fake_post)

    assert utils.get_health()["status"] == "ok"
    utils.annotate_cluster({"cluster": 1})
    utils.annotate_batch({"clusters": []})

    assert calls["get"][0].endswith("/health")
    assert calls["post"][0][0].endswith("/annotate_cluster")
    assert calls["post"][1][0].endswith("/annotate_batch")
