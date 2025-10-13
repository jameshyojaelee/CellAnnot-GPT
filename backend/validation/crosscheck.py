"""Validation utilities that cross-reference LLM annotations with the marker DB."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd


@dataclass
class CrosscheckResult:
    """Outcome of validating an LLM-proposed annotation."""

    cluster_id: str
    primary_label: str
    ontology_id: Optional[str]
    is_supported: bool
    supporting_markers: List[str] = field(default_factory=list)
    missing_markers: List[str] = field(default_factory=list)
    contradictory_markers: Dict[str, List[str]] = field(default_factory=dict)
    ontology_mismatch: bool = False
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
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
        }


def _normalize_markers(markers: Iterable[str]) -> List[str]:
    return sorted({m.upper().strip() for m in markers if isinstance(m, str) and m.strip()})


def crosscheck_annotation(
    annotation: Dict[str, Any],
    marker_db: pd.DataFrame,
    *,
    species: Optional[str] = None,
    tissue: Optional[str] = None,
    min_support: int = 1,
) -> CrosscheckResult:
    """Compare a single annotation against the reference marker database."""

    cluster_id = str(annotation.get("cluster_id", "unknown"))
    primary_label = annotation.get("primary_label", "")
    ontology_id = annotation.get("ontology_id") or None
    if not primary_label:
        raise ValueError("Annotation missing primary_label")

    markers = _normalize_markers(annotation.get("markers") or [])

    cell_types = marker_db["cell_type"].fillna("").str.lower()
    label_mask = cell_types == primary_label.lower()
    if species:
        label_mask &= marker_db["species"].fillna("").str.lower() == species.lower()
    if tissue:
        label_mask &= marker_db["tissue"].fillna("").str.lower() == tissue.lower()

    target_df = marker_db[label_mask]

    supporting: List[str] = []
    contradictory: Dict[str, List[str]] = {}

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

    notes: List[str] = []
    if not target_df.empty and not supporting:
        notes.append("Label present in DB but markers show no overlap")
    if not markers:
        notes.append("No markers supplied for validation")
    if target_df.empty:
        notes.append("Label absent from marker database")

    is_supported = (
        len(supporting) >= min_support and not contradictory and not ontology_mismatch
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
    )


def crosscheck_batch(
    annotations: Iterable[Dict[str, Any]],
    marker_db: pd.DataFrame,
    *,
    species: Optional[str] = None,
    tissue: Optional[str] = None,
    min_support: int = 1,
) -> Dict[str, CrosscheckResult]:
    """Validate a sequence of annotations and return per-cluster result mapping."""

    results: Dict[str, CrosscheckResult] = {}
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
