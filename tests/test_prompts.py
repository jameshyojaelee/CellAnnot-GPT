from __future__ import annotations

from backend.llm import prompts


def test_single_cluster_prompt_includes_candidates():
    cluster = {
        "cluster_id": "7",
        "markers": ["MS4A1"],
        "retrieved_candidates": [
            {
                "cell_type": "B cell",
                "ontology_id": "CL:0000236",
                "overlap": 2,
                "supporting_markers": ["MS4A1", "CD79A"],
                "tissue_counts": {"blood": 2},
            }
        ],
    }
    prompt = prompts.build_single_cluster_prompt(cluster, {"species": "Homo sapiens"})
    assert "Knowledge base candidates" in prompt
    assert "B cell" in prompt


def test_batch_prompt_formats_candidates():
    clusters = [
        {
            "cluster_id": "0",
            "markers": ["CD3D"],
            "retrieved_candidates": [],
        }
    ]
    prompt = prompts.build_batch_prompt(clusters, {"species": "Homo sapiens"})
    assert "Knowledge base candidates" in prompt
