"""Gene normalisation utilities supporting synonyms and ortholog mapping."""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from functools import lru_cache
from pathlib import Path
from typing import Any

from config.settings import get_settings


class GeneNormalizer:
    """Resolve gene aliases and ortholog mappings towards a primary species (default: human)."""

    def __init__(
        self,
        synonyms_path: str | Path,
        *,
        enable_orthologs: bool = True,
        primary_species: str = "Homo sapiens",
        ortholog_path: str | Path | None = None,
    ) -> None:
        self.synonyms_path = Path(synonyms_path)
        self.primary_species = primary_species
        self.primary_species_lower = primary_species.lower()
        self.enable_orthologs = enable_orthologs
        self.ortholog_path = Path(ortholog_path) if ortholog_path else None

        self.synonyms_by_species = self._load_synonyms()
        self.primary_synonyms = self.synonyms_by_species.get(self.primary_species, {})
        self.ortholog_map = self._load_orthologs() if enable_orthologs else {}

    # --------------------------------------------------------------------- load

    def _load_synonyms(self) -> dict[str, dict[str, list[str]]]:
        if not self.synonyms_path.exists():
            return {}
        with self.synonyms_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
        # Ensure consistent casing of species keys
        normalised: dict[str, dict[str, list[str]]] = {}
        for species, genes in data.items():
            species_key = str(species)
            normalised[species_key] = {gene: aliases for gene, aliases in genes.items()}
        return normalised

    def _load_orthologs(self) -> dict[str, dict[str, str]]:
        mapping: dict[str, dict[str, str]] = {}
        if not self.ortholog_path or not self.ortholog_path.exists():
            return mapping

        current_species: str | None = None
        with self.ortholog_path.open("r", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if not row:
                    continue
                if row[0].startswith("#"):
                    header = row[0].strip().lstrip("#").strip()
                    if header.lower().startswith("species="):
                        current_species = header.split("=", 1)[1].strip()
                    continue
                if len(row) < 2:
                    continue
                source_symbol, target_symbol = row[0].strip(), row[1].strip()
                if not source_symbol or not target_symbol:
                    continue
                species_key = current_species or "unknown"
                species_map = mapping.setdefault(species_key, {})
                species_map[source_symbol.upper()] = target_symbol.upper()
        return mapping

    # ------------------------------------------------------------------ helpers

    def _expand_synonyms(self, gene_upper: str, species: str) -> set[str]:
        synonyms = self.synonyms_by_species.get(species, {})
        results: set[str] = set()
        for canonical, alias_list in synonyms.items():
            canonical_upper = canonical.upper()
            alias_set = {alias.upper() for alias in alias_list}
            if gene_upper == canonical_upper or gene_upper in alias_set:
                results.add(canonical_upper)
                results.update(alias_set)
        return results

    def _canonical_in_species(self, gene_upper: str, species: str) -> str | None:
        synonyms = self.synonyms_by_species.get(species, {})
        for canonical, alias_list in synonyms.items():
            canonical_upper = canonical.upper()
            if gene_upper == canonical_upper:
                return canonical_upper
            if gene_upper in (alias.upper() for alias in alias_list):
                return canonical_upper
        if species.lower() == self.primary_species_lower:
            return gene_upper
        return None

    def _map_marker(self, marker: str, species: str) -> tuple[str, str]:
        """Return (mapped_marker_upper, status)."""

        if not marker or not isinstance(marker, str):
            return "", "unmapped"
        gene_upper = marker.upper()

        if species.lower() == self.primary_species_lower:
            canonical = self._canonical_in_species(gene_upper, self.primary_species) or gene_upper
            status = "exact" if canonical == gene_upper else "synonym"
            return canonical, status

        # Attempt species-specific synonyms first
        species_canonical = self._canonical_in_species(gene_upper, species)
        if species_canonical is None:
            species_canonical = gene_upper

        # Ortholog conversion
        if self.enable_orthologs:
            species_map = self.ortholog_map.get(species, {})
            mapped = species_map.get(species_canonical) or species_map.get(gene_upper)
            if mapped:
                return mapped, "ortholog"

        # If already human-like, return canonical human form
        human_canonical = self._canonical_in_species(gene_upper, self.primary_species)
        if human_canonical:
            status = "presumed_human" if human_canonical == gene_upper else "synonym"
            return human_canonical, status

        return gene_upper, "unmapped"

    # ---------------------------------------------------------- public methods

    def normalise_marker(self, gene: str, species: str | None = None) -> set[str]:
        """Return the set of recognised aliases for a gene within a species."""

        if not gene or not isinstance(gene, str):
            return set()

        species_name = species or self.primary_species
        gene_upper = gene.upper()
        results: set[str] = {gene_upper}
        results.update(self._expand_synonyms(gene_upper, species_name))

        if species_name.lower() != self.primary_species_lower:
            # Also include primary-species synonyms if ortholog mapping lands there
            human_synonyms = self._expand_synonyms(gene_upper, self.primary_species)
            results.update(human_synonyms)

        return results

    def normalise_markers(self, markers: Iterable[str], species: str | None = None) -> list[str]:
        """Return a deduplicated list of normalised markers for a given species."""

        seen: set[str] = set()
        ordered: list[str] = []
        for marker in markers:
            for alias in self.normalise_marker(marker, species):
                if alias and alias not in seen:
                    seen.add(alias)
                    ordered.append(alias)
        return ordered

    def map_to_primary(
        self,
        markers: Iterable[str],
        species: str | None = None,
    ) -> tuple[list[str], list[dict[str, Any]]]:
        """Map markers to the primary species, returning mapped markers and mapping notes."""

        species_name = species or self.primary_species
        mapped_markers: list[str] = []
        notes: list[dict[str, Any]] = []

        for marker in markers:
            mapped, status = self._map_marker(marker, species_name)
            if not mapped:
                continue
            mapped_markers.append(mapped)
            if status != "exact":
                notes.append(
                    {
                        "source": marker,
                        "target": mapped if status != "unmapped" else None,
                        "status": status,
                        "species": species_name,
                    }
                )

        if not mapped_markers:
            mapped_markers = [
                marker.upper() for marker in markers if isinstance(marker, str) and marker.strip()
            ]

        return mapped_markers, notes


@lru_cache(maxsize=1)
def get_gene_normalizer() -> GeneNormalizer:
    settings = get_settings()
    return GeneNormalizer(
        settings.synonym_config_path,
        enable_orthologs=settings.synonym_enable_orthologs,
        primary_species=settings.kb_primary_species,
        ortholog_path=settings.ortholog_mapping_path,
    )


__all__ = ["GeneNormalizer", "get_gene_normalizer"]
