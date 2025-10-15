"""Prompt builders for GPT Cell Annotator's LLM interactions."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any


def build_system_prompt() -> str:
    """Return the default system prompt that constrains model behaviour."""

    return (
        "You are GPT Cell Annotator, an expert single-cell biologist and computational analyst. "
        "Always ground annotations in marker-gene evidence, cite databases such as PanglaoDB, "
        "CellMarker, or relevant literature, and surface uncertainty transparently. "
        "Prefer clear JSON outputs matching the requested schema."
    )


def _format_context(dataset_context: dict[str, Any] | None) -> str:
    if not dataset_context:
        return "Context: not provided."
    parts = [f"{key}: {value}" for key, value in dataset_context.items() if value]
    return "Context\\n- " + "\\n- ".join(parts) if parts else "Context: not provided."


def _format_candidates(cluster_payload: dict[str, Any]) -> str:
    candidates = cluster_payload.get("retrieved_candidates") or []
    if not candidates:
        return "Knowledge base candidates: none retrieved."
    lines: list[str] = []
    for idx, candidate in enumerate(candidates, start=1):
        supporting = ", ".join(candidate.get("supporting_markers") or []) or "n/a"
        ontology_id = candidate.get("ontology_id") or "unknown"
        overlap = candidate.get("overlap", 0)
        lines.append(
            f"{idx}. {candidate.get('cell_type', 'Unknown')} (ontology_id={ontology_id}, overlap={overlap}) supporting markers: {supporting}"
        )
    return "Knowledge base candidates:\\n" + "\\n".join(f"- {line}" for line in lines)


def build_single_cluster_prompt(
    cluster_payload: dict[str, Any], dataset_context: dict[str, Any] | None = None
) -> str:
    """Construct a single-cluster annotation prompt."""

    context_block = _format_context(dataset_context)
    markers = cluster_payload.get("markers", [])
    marker_lines = "\\n".join(f"- {gene}" for gene in markers) or "- No markers supplied."
    candidate_block = _format_candidates(cluster_payload)

    stats = cluster_payload.get("stats") or {}
    stats_lines = "\\n".join(f"{key}: {value}" for key, value in stats.items())
    stats_block = (
        f"Cluster stats:\\n{stats_lines}" if stats_lines else "Cluster stats: not provided."
    )

    cluster_id = cluster_payload.get("cluster_id", "unknown")

    instructions = (
        "Return JSON with keys: primary_label, ontology_id, confidence, rationale, "
        "alternatives, caveats, evidence_sources. Provide alternatives as a list with 'label' "
        "and 'reason'. Use 'Unknown or Novel' if evidence is conflicting. "
        "Evidence_sources must list at least one citation referencing PanglaoDB, CellMarker, "
        "or curated literature. Explicitly evaluate the knowledge-base candidates above, noting "
        "why each matches or conflicts and copy ontology IDs from the evidence when available."
    )

    return (
        f"{context_block}\\n"
        f"Cluster ID: {cluster_id}\\n"
        f"Marker genes (descending importance):\\n{marker_lines}\\n"
        f"{candidate_block}\\n"
        f"{stats_block}\\n"
        f"{instructions}"
    )


def build_batch_prompt(
    clusters: Iterable[dict[str, Any]],
    dataset_context: dict[str, Any] | None = None,
) -> str:
    """Construct a batch annotation prompt for multiple clusters."""

    context_block = _format_context(dataset_context)
    cluster_sections = []
    for cluster in clusters:
        cluster_id = cluster.get("cluster_id", "unknown")
        markers = cluster.get("markers", [])
        marker_lines = "\\n    ".join(markers) or "    No markers supplied."
        stats = cluster.get("stats") or {}
        stats_lines = ", ".join(f"{k}: {v}" for k, v in stats.items()) or "None"
        candidates_block = _format_candidates(cluster).replace("\\n", "\\n    ")
        cluster_sections.append(
            f"- Cluster {cluster_id}:\\n"
            f"    markers:\\n    {marker_lines}\\n"
            f"    stats: {stats_lines}\\n"
            f"    {candidates_block}\\n"
            f"    notes: {cluster.get('notes', 'None')}"
        )

    cluster_block = "\\n".join(cluster_sections)
    instructions = (
        "For each cluster, return an object with primary_label, ontology_id, confidence, "
        "rationale, alternatives, caveats, evidence_sources (with at least one knowledge-base "
        "citation such as PanglaoDB, CellMarker, or curated literature). Explicitly evaluate the "
        "knowledge-base candidates, stating why they support, conflict with, or fail to explain "
        "the cluster. Provide an overall dataset_summary with observed composition, uncertainties, "
        "and follow-up recommendations."
    )

    return f"{context_block}\\nClusters:\\n{cluster_block}\\n{instructions}"


def build_uncertainty_prompt(
    cluster_payload: dict[str, Any], dataset_context: dict[str, Any] | None = None
) -> str:
    """Construct a prompt that analyses why annotation is ambiguous."""

    context_block = _format_context(dataset_context)
    markers = cluster_payload.get("markers", [])
    marker_lines = "\\n".join(f"- {gene}" for gene in markers) or "- No markers supplied."
    current_guess = cluster_payload.get("current_guess", "Unknown")

    instructions = (
        "1. Enumerate competing hypotheses with pros/cons based on markers and tissue.\\n"
        "2. Suggest follow-up experiments or markers to resolve ambiguity.\\n"
        "3. Provide a cautious report sentence acknowledging uncertainty."
    )

    return (
        f"{context_block}\\n"
        f"Markers:\\n{marker_lines}\\n"
        f"Existing model guess: {current_guess}\\n"
        f"{instructions}"
    )


__all__ = [
    "build_system_prompt",
    "build_single_cluster_prompt",
    "build_batch_prompt",
    "build_uncertainty_prompt",
]

