from __future__ import annotations

from backend.util.gene_normalization import GeneNormalizer, get_gene_normalizer


def test_gene_normaliser_synonyms(tmp_path) -> None:
    config_path = tmp_path / "synonyms.json"
    config_path.write_text(
        """
        {
          "Homo sapiens": {
            "MS4A1": ["CD20"]
          },
          "orthologs": {
            "Ms4a1": "MS4A1"
          }
        }
        """,
        encoding="utf-8",
    )
    normalizer = GeneNormalizer(config_path, enable_orthologs=True)
    assert {"MS4A1", "CD20"} <= normalizer.normalise_marker("MS4A1", "Homo sapiens")
    assert "MS4A1" in normalizer.normalise_marker("CD20", "Homo sapiens")
    mouse = normalizer.normalise_marker("Ms4a1", "Mus musculus")
    assert "MS4A1" in mouse


def test_global_gene_normaliser_uses_settings() -> None:
    normalizer = get_gene_normalizer()
    assert "MS4A1" in normalizer.normalise_marker("cd20", "Homo sapiens")
