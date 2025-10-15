from __future__ import annotations

import json
from typing import Any

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from backend.api import main
from backend.api.main import app


class DummyAnnotator:
    def __init__(self) -> None:
        self.cluster_calls = 0
        self.batch_calls = 0
        self._llm_mode = "mock"

    def annotate_cluster(self, cluster_payload, dataset_context=None):  # type: ignore[override]
        self.cluster_calls += 1
        return {
            "primary_label": "B cell",
            "ontology_id": "CL:0000236",
            "confidence": "High",
            "rationale": "MS4A1 matches B cell markers.",
        }

    def annotate_batch(self, clusters, dataset_context=None):  # type: ignore[override]
        self.batch_calls += 1
        return {
            str(c["cluster_id"]): {
                "primary_label": "T cell",
                "ontology_id": "CL:0000084",
                "confidence": "Medium",
                "rationale": "CD3E is a canonical T cell marker.",
            }
            for c in clusters
        }

    @property
    def llm_mode(self) -> str:  # type: ignore[override]
        return self._llm_mode


@pytest.fixture(autouse=True)
def override_dependencies():
    annotator = DummyAnnotator()
    marker_db = pd.DataFrame(
        [
            {
                "source": "PanglaoDB",
                "cell_type": "B cell",
                "ontology_id": "CL:0000236",
                "gene_symbol": "MS4A1",
                "species": "Homo sapiens",
                "tissue": "Blood",
                "evidence": "",
                "reference": "",
            },
            {
                "source": "PanglaoDB",
                "cell_type": "T cell",
                "ontology_id": "CL:0000084",
                "gene_symbol": "CD3E",
                "species": "Homo sapiens",
                "tissue": "Blood",
                "evidence": "",
                "reference": "",
            },
        ]
    )

    main.app.dependency_overrides[main.get_annotator] = lambda: annotator
    main.app.dependency_overrides[main.get_marker_db] = lambda: marker_db
    main.app.dependency_overrides[main.get_cache] = lambda: None
    yield
    main.app.dependency_overrides.clear()


client = TestClient(app)


class FakeCache:
    def __init__(self) -> None:
        self.store: dict[str, Any] = {}
        self.get_calls = 0
        self.set_calls = 0

    def _key(self, payload: dict[str, Any]) -> str:
        return json.dumps(payload, sort_keys=True)

    async def get(self, payload: dict[str, Any]) -> dict[str, Any] | None:
        self.get_calls += 1
        return self.store.get(self._key(payload))

    async def set(self, payload: dict[str, Any], value: dict[str, Any]) -> None:
        self.set_calls += 1
        self.store[self._key(payload)] = value


class FailingCache:
    async def get(self, payload: dict[str, Any]) -> None:
        raise RuntimeError("redis unavailable")

    async def set(self, payload: dict[str, Any], value: dict[str, Any]) -> None:
        raise RuntimeError("redis unavailable")


def test_health_endpoint():
    response = client.get("/health")
    assert response.status_code == 200
    assert response.json() == {"status": "ok", "llm_mode": "mock", "cache_enabled": False}


def test_annotate_cluster_endpoint():
    payload = {
        "cluster": {
            "cluster_id": "0",
            "markers": ["MS4A1"],
        },
        "dataset_context": {"species": "Homo sapiens", "tissue": "Blood"},
        "return_validated": True,
    }
    response = client.post("/annotate_cluster", json=payload)
    assert response.status_code == 200
    data = response.json()
    assert data["result"]["summary"]["supported_clusters"] == 1
    assert "support_rate" in data["result"]["metrics"]
    assert data["result"]["clusters"][0]["annotation"]["primary_label"] == "B cell"
    assert data["result"]["clusters"][0]["status"] == "supported"


def test_annotate_batch_endpoint():
    payload = {
        "clusters": [
            {"cluster_id": "0", "markers": ["CD3E"]},
            {"cluster_id": "1", "markers": ["CD3E"]},
        ],
        "dataset_context": {"species": "Homo sapiens"},
    }
    response = client.post("/annotate_batch", json=payload)
    assert response.status_code == 200
    data = response.json()
    assert data["result"]["summary"]["total_clusters"] == 2
    assert {cluster["annotation"]["primary_label"] for cluster in data["result"]["clusters"]} == {
        "T cell"
    }
    assert data["result"]["metrics"]["confidence_counts"]


def test_annotate_cluster_without_validation_returns_raw():
    payload = {
        "cluster": {"cluster_id": "0", "markers": ["MS4A1"]},
        "dataset_context": {"species": "Homo sapiens"},
        "return_validated": False,
    }
    response = client.post("/annotate_cluster", json=payload)
    assert response.status_code == 200
    result = response.json()["result"]
    assert result["primary_label"] == "B cell"
    assert result["cluster_id"] == "0"
    assert "summary" not in result


def test_annotate_cluster_uses_cache_for_non_validated():
    fake_cache = FakeCache()
    main.app.dependency_overrides[main.get_cache] = lambda: fake_cache
    payload = {
        "cluster": {"cluster_id": "0", "markers": ["MS4A1"]},
        "dataset_context": {"species": "Homo sapiens"},
        "return_validated": False,
    }
    first = client.post("/annotate_cluster", json=payload)
    assert first.status_code == 200
    assert fake_cache.set_calls == 1
    assert fake_cache.get_calls == 1
    cache_key = next(iter(fake_cache.store.keys()))
    assert '"validated": false' in cache_key
    stored_payload = next(iter(fake_cache.store.values()))
    assert stored_payload["primary_label"] == "B cell"
    assert stored_payload["cluster_id"] == "0"

    second = client.post("/annotate_cluster", json=payload)
    assert second.status_code == 200
    assert fake_cache.get_calls == 2
    assert fake_cache.set_calls == 1  # second call served from cache
    assert second.json()["result"] == first.json()["result"]


def test_annotate_cluster_validated_bypasses_unvalidated_cache():
    fake_cache = FakeCache()
    main.app.dependency_overrides[main.get_cache] = lambda: fake_cache

    # prime cache with non-validated response
    non_validated_payload = {
        "cluster": {"cluster_id": "0", "markers": ["MS4A1"]},
        "dataset_context": {"species": "Homo sapiens"},
        "return_validated": False,
    }
    client.post("/annotate_cluster", json=non_validated_payload)
    assert fake_cache.set_calls == 1

    validated_payload = {
        "cluster": {"cluster_id": "0", "markers": ["MS4A1"]},
        "dataset_context": {"species": "Homo sapiens"},
        "return_validated": True,
    }
    response = client.post("/annotate_cluster", json=validated_payload)
    assert response.status_code == 200
    # validated path should not reuse cached raw entry
    assert response.json()["result"]["summary"]["total_clusters"] == 1


def test_cache_failures_do_not_break_endpoint():
    main.app.dependency_overrides[main.get_cache] = lambda: FailingCache()
    payload = {
        "cluster": {"cluster_id": "0", "markers": ["MS4A1"]},
        "return_validated": False,
    }
    response = client.post("/annotate_cluster", json=payload)
    assert response.status_code == 200
    assert response.json()["result"]["primary_label"] == "B cell"
