from __future__ import annotations

import pandas as pd

from backend.validation.crosscheck import (
    CrosscheckResult,
    crosscheck_annotation,
    crosscheck_batch,
)
from backend.validation.report import build_structured_report, render_text_report


def make_marker_db() -> pd.DataFrame:
    return pd.DataFrame(
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
            {
                "source": "CellMarker",
                "cell_type": "NK cell",
                "ontology_id": "CL:0000623",
                "gene_symbol": "GNLY",
                "species": "Homo sapiens",
                "tissue": "Blood",
                "evidence": "",
                "reference": "",
            },
        ]
    )


def test_crosscheck_detects_contradiction() -> None:
    annotation = {
        "cluster_id": "0",
        "primary_label": "B cell",
        "ontology_id": "CL:0000236",
        "markers": ["MS4A1", "GNLY"],
    }
    result = crosscheck_annotation(
        annotation,
        make_marker_db(),
        species="Homo sapiens",
        tissue="Blood",
    )

    assert result.is_supported is False
    assert result.supporting_markers == ["MS4A1"]
    assert "GNLY" in result.contradictory_markers
    assert result.ontology_mismatch is False


def test_crosscheck_supports_valid_annotation() -> None:
    annotation = {
        "cluster_id": "1",
        "primary_label": "T cell",
        "ontology_id": "CL:0000084",
        "markers": ["CD3E"],
    }
    result = crosscheck_annotation(
        annotation,
        make_marker_db(),
        species="Homo sapiens",
    )

    assert result.is_supported is True
    assert result.supporting_markers == ["CD3E"]
    assert result.contradictory_markers == {}
    assert result.missing_markers == []


def test_crosscheck_handles_unknown_label() -> None:
    annotation = {
        "cluster_id": "2",
        "primary_label": "Monocyte",
        "ontology_id": "CL:0000576",
        "markers": ["LYZ"],
    }
    result = crosscheck_annotation(annotation, make_marker_db())

    assert result.is_supported is False
    assert "Label absent from marker database" in result.notes


def test_crosscheck_batch_returns_mapping() -> None:
    annotations = [
        {"cluster_id": "0", "primary_label": "B cell", "markers": ["MS4A1"]},
        {"cluster_id": "1", "primary_label": "T cell", "markers": ["CD3E"]},
    ]
    results = crosscheck_batch(annotations, make_marker_db(), species="Homo sapiens")
    assert isinstance(results["0"], CrosscheckResult)
    assert results["0"].primary_label == "B cell"
    assert results["1"].is_supported is True


def test_build_structured_report_and_render_text() -> None:
    annotations = [
        {
            "cluster_id": "0",
            "primary_label": "B cell",
            "confidence": "High",
            "rationale": "Classic B cell markers",
            "markers": ["MS4A1"],
        },
        {
            "cluster_id": "1",
            "primary_label": "Unknown or Novel",
            "confidence": "Low",
            "rationale": "Markers do not match known profiles",
            "markers": ["GNLY"],
        },
    ]
    results = crosscheck_batch(annotations, make_marker_db(), species="Homo sapiens")
    structured = build_structured_report(annotations, results)

    assert structured.summary.total_clusters == 2
    assert structured.summary.flagged_clusters >= 1
    assert "1" in structured.summary.unknown_clusters
    assert "confidence_counts" in structured.metrics.model_dump()

    text = render_text_report(structured)
    assert "GPT Cell Annotator Validation Report" in text
    assert "Unknown clusters" in text
