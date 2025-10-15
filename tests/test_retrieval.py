from __future__ import annotations

import sqlite3
from pathlib import Path

from backend.llm import retrieval
from backend.llm.retrieval import MarkerRetriever, retrieve_candidates
from config.settings import Settings, get_settings


def _build_sqlite(db_path: Path) -> None:
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute(
        """
        CREATE TABLE cell_markers (
            source TEXT,
            cell_type TEXT,
            ontology_id TEXT,
            gene_symbol TEXT,
            species TEXT,
            tissue TEXT,
            evidence TEXT,
            reference TEXT
        )
        """
    )
    cur.executemany(
        """
        INSERT INTO cell_markers (source, cell_type, ontology_id, gene_symbol, species, tissue, evidence, reference)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            ("Demo", "B cell", "CL:0000236", "MS4A1", "Homo sapiens", "Blood", "", ""),
            ("Demo", "B cell", "CL:0000236", "CD79A", "Homo sapiens", "Blood", "", ""),
            ("Demo", "T cell", "CL:0000084", "CD3D", "Homo sapiens", "Blood", "", ""),
            ("Demo", "T cell", "CL:0000084", "CD3E", "Homo sapiens", "Blood", "", ""),
        ],
    )
    conn.commit()
    conn.close()


def test_marker_retriever_returns_candidates(tmp_path: Path) -> None:
    db_path = tmp_path / "marker_db.sqlite"
    _build_sqlite(db_path)
    retriever = MarkerRetriever(db_path)

    candidates = retriever.retrieve(
        ["MS4A1", "CD79A"],
        species="Homo sapiens",
        tissue="Blood",
        top_k=5,
        min_overlap=1,
    )
    assert candidates
    assert candidates[0].cell_type == "B cell"
    assert "MS4A1" in candidates[0].supporting_markers


def test_retrieve_candidates_integration(tmp_path: Path, monkeypatch) -> None:
    db_path = tmp_path / "marker_db.sqlite"
    _build_sqlite(db_path)

    current = Settings()

    class DummySettings(Settings):
        def __init__(self) -> None:
            super().__init__()
            self.data_dir = str(tmp_path)
            self.rag_enabled = True
            self.rag_top_k = 3
            self.rag_min_overlap = 1

    dummy = DummySettings()
    monkeypatch.setattr("backend.llm.retrieval.get_settings", lambda: dummy)
    retrieval.get_retriever.cache_clear()  # type: ignore[attr-defined]

    candidates = retrieve_candidates(
        ["CD3d"], species="Homo sapiens", tissue="Blood"
    )
    assert candidates
    assert candidates[0].cell_type == "T cell"
