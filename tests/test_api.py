from __future__ import annotations

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
