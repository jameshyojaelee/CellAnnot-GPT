from __future__ import annotations

import pandas as pd
import pytest

from backend.validation.crosscheck import crosscheck_annotation
from backend.validation.report import build_structured_report
from config.settings import get_settings
from frontend.streamlit_app import _parse_marker_field


def make_marker_db() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "source": "Demo",
                "cell_type": "B cell",
                "ontology_id": "CL:0000236",
                "gene_symbol": "MS4A1",
                "species": "Homo sapiens",
                "tissue": "Blood",
                "evidence": "",
                "reference": "",
                "evidence_score": "high",
            }
        ]
    )


def test_build_report_for_single_cluster():
    get_settings().validation_min_marker_overlap = 1
    annotation = {
        "cluster_id": "single",
        "primary_label": "B cell",
        "ontology_id": "CL:0000236",
        "confidence": "High",
        "rationale": "MS4A1 matches B cell markers",
        "markers": ["MS4A1"],
    }
    result = crosscheck_annotation(annotation, make_marker_db(), species="Homo sapiens")
    report = build_structured_report([annotation], {"single": result})

    assert report.summary.total_clusters == 1
    assert report.summary.supported_clusters == 1
    assert report.metrics.support_rate == 1.0
    assert report.clusters[0].status == "supported"


@pytest.mark.parametrize(
    "value,expected",
    [
        ("MS4A1, CD79A", ["MS4A1", "CD79A"]),
        ('["MS4A1", "CD79A"]', ["MS4A1", "CD79A"]),
        (["MS4A1", "CD79A"], ["MS4A1", "CD79A"]),
    ],
)
def test_parse_marker_field(value, expected):
    assert _parse_marker_field(value) == expected
