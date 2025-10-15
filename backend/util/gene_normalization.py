"""Gene normalisation utilities supporting synonym and ortholog mapping."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Set

from config.settings import get_settings


class GeneNormalizer:
    """Resolve gene aliases and optional ortholog mappings."""

    def __init__(self, config_path: str | Path, enable_orthologs: bool = True) -> None:
        self.config_path = Path(config_path)
        self.enable_orthologs = enable_orthologs
        self._synonyms = self._load_synonyms()

    def _load_synonyms(self) -> dict[str, dict[str, list[str]]]:
        if not self.config_path.exists():
            return {}
        with self.config_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
        return data

    @lru_cache(maxsize=512)
    def normalise_marker(self, gene: str, species: str | None = None) -> Set[str]:
        """Return a set of canonical markers for the provided gene."""

        if not gene or not isinstance(gene, str):
            return set()

        upper_gene = gene.upper()
        results: set[str] = {upper_gene}
        human_synonyms = self._synonyms.get("Homo sapiens", {})

        if species:
            species_synonyms = self._synonyms.get(species, {})
            if species_synonyms:
                results.update(alias.upper() for alias in species_synonyms.get(gene, []))
                results.update(alias.upper() for alias in species_synonyms.get(upper_gene, []))
                for canonical, alias_list in species_synonyms.items():
                    if upper_gene == canonical.upper() or upper_gene in {alias.upper() for alias in alias_list}:
                        results.add(canonical.upper())
                        results.update(alias.upper() for alias in alias_list)

        human_aliases = human_synonyms.get(upper_gene, [])
        results.update(alias.upper() for alias in human_aliases)
        for canonical, alias_list in human_synonyms.items():
            alias_set = {alias.upper() for alias in alias_list}
            if upper_gene == canonical.upper() or upper_gene in alias_set:
                results.add(canonical.upper())
                results.update(alias_set)

        if self.enable_orthologs:
            orthologs = self._synonyms.get("orthologs", {})
            human_equivalent = orthologs.get(gene) or orthologs.get(upper_gene)
            if human_equivalent:
                results.add(human_equivalent.upper())
                results.update(alias.upper() for alias in human_synonyms.get(human_equivalent.upper(), []))

        return results

    def normalise_markers(self, markers: Iterable[str], species: str | None = None) -> Set[str]:
        results: set[str] = set()
        for marker in markers:
            results.update(self.normalise_marker(marker, species))
        return results


@lru_cache(maxsize=1)
def get_gene_normalizer() -> GeneNormalizer:
    settings = get_settings()
    return GeneNormalizer(
        settings.synonym_config_path,
        enable_orthologs=settings.synonym_enable_orthologs,
    )


__all__ = ["GeneNormalizer", "get_gene_normalizer"]
