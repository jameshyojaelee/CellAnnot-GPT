"""Knowledge-base retrieval helpers for prompt augmentation."""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

from backend.util.gene_normalization import get_gene_normalizer
from config.settings import get_settings
from gpt_cell_annotator import assets


@dataclass(frozen=True)
class RetrievalCandidate:
    cell_type: str
    ontology_id: str | None
    overlap: int
    supporting_markers: list[str]
    tissue_counts: dict[str, int]


class MarkerRetriever:
    """Query the marker sqlite database for overlapping cell types."""

    def __init__(self, sqlite_path: str | Path) -> None:
        self.sqlite_path = Path(sqlite_path)
        self._connection: sqlite3.Connection | None = None

    def _connect(self) -> sqlite3.Connection:
        if self._connection is None:
            if not self.sqlite_path.exists():
                raise FileNotFoundError(f"Marker sqlite DB not found at {self.sqlite_path}")
            self._connection = sqlite3.connect(str(self.sqlite_path))
        return self._connection

    def close(self) -> None:
        if self._connection is not None:
            self._connection.close()
            self._connection = None

    def retrieve(
        self,
        markers: Iterable[str],
        *,
        species: str | None,
        tissue: str | None,
        top_k: int,
        min_overlap: int,
    ) -> list[RetrievalCandidate]:
        connection = self._connect()
        normalizer = get_gene_normalizer()
        expanded_markers = normalizer.normalise_markers(markers, species)
        if not expanded_markers:
            return []

        placeholders = ",".join("?" for _ in expanded_markers)
        params: list[str] = list(expanded_markers)
        species_filter = ""
        tissue_filter = ""
        if species:
            species_filter = "AND LOWER(species) = LOWER(?)"
            params.append(species)
        if tissue:
            tissue_filter = "AND LOWER(tissue) = LOWER(?)"
            params.append(tissue)

        query = f"""
            SELECT cell_type, ontology_id, gene_symbol, tissue
            FROM cell_markers
            WHERE UPPER(gene_symbol) IN ({placeholders})
            {species_filter}
            {tissue_filter}
        """

        cursor = connection.execute(query, params)
        rows = cursor.fetchall()
        if not rows:
            return []

        counts: dict[str, dict[str, object]] = {}
        for cell_type, ontology_id, gene_symbol, tissue_value in rows:
            cell_key = cell_type or "Unknown"
            entry = counts.setdefault(
                cell_key,
                {
                    "ontology_id": ontology_id or None,
                    "genes": set(),
                    "tissue_counts": {},
                },
            )
            entry["genes"].add(str(gene_symbol).upper())
            tissues = entry["tissue_counts"]
            tissue_key = (tissue_value or "unspecified").lower()
            tissues[tissue_key] = tissues.get(tissue_key, 0) + 1

        candidates: list[RetrievalCandidate] = []
        for cell_type, info in counts.items():
            supporting_genes = sorted(info["genes"])
            overlap = len(supporting_genes)
            if overlap < min_overlap:
                continue
            candidates.append(
                RetrievalCandidate(
                    cell_type=cell_type,
                    ontology_id=info["ontology_id"],
                    overlap=overlap,
                    supporting_markers=supporting_genes,
                    tissue_counts=dict(info["tissue_counts"]),
                )
            )

        candidates.sort(key=lambda item: (item.overlap, item.cell_type), reverse=True)
        return candidates[:top_k]


@lru_cache(maxsize=1)
def get_retriever() -> MarkerRetriever | None:
    settings = get_settings()
    data_dir = Path(settings.data_dir)
    sqlite_path = data_dir / "marker_db.sqlite"

    if not sqlite_path.exists():
        assets.ensure_marker_database(target_dir=data_dir)
        sqlite_path = data_dir / "marker_db.sqlite"

    if sqlite_path.exists():
        return MarkerRetriever(sqlite_path)

    try:
        packaged_dir = assets.ensure_marker_database()
    except FileNotFoundError:
        return None
    fallback_path = packaged_dir / "marker_db.sqlite"
    if not fallback_path.exists():
        return None
    data_dir.mkdir(parents=True, exist_ok=True)
    sqlite_target = data_dir / "marker_db.sqlite"
    sqlite_target.write_bytes(fallback_path.read_bytes())
    return MarkerRetriever(sqlite_target)


def retrieve_candidates(
    markers: Sequence[str],
    *,
    species: str | None,
    tissue: str | None,
) -> list[RetrievalCandidate]:
    settings = get_settings()
    if not settings.rag_enabled:
        return []
    retriever = get_retriever()
    if retriever is None:
        return []
    try:
        return retriever.retrieve(
            markers,
            species=species,
            tissue=tissue,
            top_k=settings.rag_top_k,
            min_overlap=settings.rag_min_overlap,
        )
    except Exception:
        return []


__all__ = ["retrieve_candidates", "RetrievalCandidate", "get_retriever", "MarkerRetriever"]
