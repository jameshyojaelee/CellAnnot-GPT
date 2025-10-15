"""Validation utilities that cross-reference LLM annotations with the marker DB."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import pandas as pd


@dataclass
class CrosscheckResult:
    """Outcome of validating an LLM-proposed annotation."""

    cluster_id: str
    primary_label: str
    ontology_id: str | None
    is_supported: bool
    supporting_markers: list[str] = field(default_factory=list)
    missing_markers: list[str] = field(default_factory=list)
    contradictory_markers: dict[str, list[str]] = field(default_factory=dict)
    ontology_mismatch: bool = False
    notes: list[str] = field(default_factory=list)
    support_count: int = 0
    flag_reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "cluster_id": self.cluster_id,
            "primary_label": self.primary_label,
            "ontology_id": self.ontology_id,
            "is_supported": self.is_supported,
            "supporting_markers": self.supporting_markers,
            "missing_markers": self.missing_markers,
            "contradictory_markers": self.contradictory_markers,
            "ontology_mismatch": self.ontology_mismatch,
            "notes": self.notes,
            "support_count": self.support_count,
            "flag_reasons": self.flag_reasons,
        }


def _normalize_markers(markers: Iterable[str]) -> list[str]:
    return sorted({m.upper().strip() for m in markers if isinstance(m, str) and m.strip()})


def crosscheck_annotation(
    annotation: dict[str, Any],
    marker_db: pd.DataFrame,
    *,
    species: str | None = None,
    tissue: str | None = None,
    min_support: int = 1,
) -> CrosscheckResult:
    """Compare a single annotation against the reference marker database."""

    cluster_id = str(annotation.get("cluster_id", "unknown"))
    primary_label = annotation.get("primary_label", "")
    ontology_id = annotation.get("ontology_id") or None
    if not primary_label:
        raise ValueError("Annotation missing primary_label")

    markers = _normalize_markers(annotation.get("markers") or [])

    # Track contextual subsets to infer reason codes.
    cell_types = marker_db["cell_type"].fillna("").str.lower()
    label_mask = cell_types == primary_label.lower()
    label_df = marker_db[label_mask]

    species_df = label_df
    species_mismatch = False
    if species:
        species_mask = species_df["species"].fillna("").str.lower() == species.lower()
        filtered = species_df[species_mask]
        if species_df.shape[0] > 0 and filtered.shape[0] == 0:
            species_mismatch = True
        species_df = filtered if not filtered.empty else species_df

    tissue_df = species_df
    tissue_mismatch = False
    if tissue:
        tissue_mask = tissue_df["tissue"].fillna("").str.lower() == tissue.lower()
        filtered = tissue_df[tissue_mask]
        if tissue_df.shape[0] > 0 and filtered.shape[0] == 0:
            tissue_mismatch = True
        tissue_df = filtered if not filtered.empty else tissue_df

    target_df = tissue_df

    flag_reasons: list[str] = []
    if label_df.empty:
        flag_reasons.append("label_not_in_kb")
    if species_mismatch:
        flag_reasons.append("species_mismatch")
    if tissue_mismatch:
        flag_reasons.append("tissue_mismatch")

    supporting: list[str] = []
    contradictory: dict[str, list[str]] = {}

    if not target_df.empty:
        target_markers = _normalize_markers(target_df["gene_symbol"].tolist())
        supporting = sorted(set(markers) & set(target_markers))

    cell_marker_map = (
        marker_db.assign(
            __gene=marker_db["gene_symbol"].fillna("").str.upper().str.strip(),
            __cell=marker_db["cell_type"].fillna("").str.strip(),
        )
        .groupby("__gene")["__cell"]
        .apply(lambda s: sorted({c for c in s if c}))
    )

    for marker in markers:
        cell_types_for_marker = cell_marker_map.get(marker, [])
        if cell_types_for_marker and primary_label not in cell_types_for_marker:
            contradictory[marker] = cell_types_for_marker

    missing = sorted(set(markers) - set(supporting) - set(contradictory.keys()))

    ontology_mismatch = False
    if ontology_id:
        matches = target_df["ontology_id"].dropna().str.upper().str.strip()
        ontology_mismatch = matches.empty or ontology_id.upper().strip() not in set(matches)

    notes: list[str] = []
    if not target_df.empty and not supporting:
        notes.append("Label present in DB but markers show no overlap")
    if not markers:
        notes.append("No markers supplied for validation")
        flag_reasons.append("no_markers_supplied")
    if target_df.empty:
        notes.append("Label absent from marker database")
    if species_mismatch:
        notes.append("No marker entries for requested species")
    if tissue_mismatch:
        notes.append("No marker entries for requested tissue")
    if not ontology_id:
        notes.append("Ontology identifier missing from annotation")
        flag_reasons.append("missing_ontology_id")

    support_count = len(supporting)
    if support_count < min_support:
        flag_reasons.append("low_marker_overlap")
    if contradictory:
        flag_reasons.append("contradictory_markers")
    if ontology_mismatch:
        flag_reasons.append("ontology_mismatch")

    is_supported = (
        support_count >= max(min_support, 1)
        and not contradictory
        and not ontology_mismatch
        and "label_not_in_kb" not in flag_reasons
        and "species_mismatch" not in flag_reasons
        and "tissue_mismatch" not in flag_reasons
        and "no_markers_supplied" not in flag_reasons
        and "missing_ontology_id" not in flag_reasons
    )

    return CrosscheckResult(
        cluster_id=cluster_id,
        primary_label=primary_label,
        ontology_id=ontology_id,
        is_supported=is_supported,
        supporting_markers=supporting,
        missing_markers=missing,
        contradictory_markers=contradictory,
        ontology_mismatch=ontology_mismatch,
        notes=notes,
        support_count=support_count,
        flag_reasons=sorted(set(flag_reasons)),
    )


def crosscheck_batch(
    annotations: Iterable[dict[str, Any]],
    marker_db: pd.DataFrame,
    *,
    species: str | None = None,
    tissue: str | None = None,
    min_support: int = 1,
) -> dict[str, CrosscheckResult]:
    """Validate a sequence of annotations and return per-cluster result mapping."""

    results: dict[str, CrosscheckResult] = {}
    for annotation in annotations:
        cluster_id = str(annotation.get("cluster_id", "unknown"))
        result = crosscheck_annotation(
            annotation,
            marker_db,
            species=species,
            tissue=tissue,
            min_support=min_support,
        )
        results[cluster_id] = result
    return results


__all__ = ["CrosscheckResult", "crosscheck_annotation", "crosscheck_batch"]
