from __future__ import annotations

import json
from pathlib import Path

import pytest

from config import settings
from gpt_cell_annotator import cli


@pytest.fixture(autouse=True)
def reset_settings_cache() -> None:
    settings.get_settings.cache_clear()


def test_cli_annotate_offline(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    assets_home = tmp_path / "assets"
    monkeypatch.setenv("GPT_CELL_ANNOTATOR_HOME", str(assets_home))
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)

    output_json = tmp_path / "annotations.json"
    exit_code = cli.main(
        [
            "annotate",
            "data/demo/pbmc_markers.csv",
            "--offline",
            "--out-json",
            str(output_json),
        ]
    )

    assert exit_code == 0
    captured = capsys.readouterr().out
    assert "Cluster" in captured
    assert output_json.exists()

    report = json.loads(output_json.read_text(encoding="utf-8"))
    assert "clusters" in report
    assert any(cluster.get("annotation", {}).get("primary_label") for cluster in report["clusters"])


def test_cli_build_db_offline(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    assets_home = tmp_path / "assets"
    monkeypatch.setenv("GPT_CELL_ANNOTATOR_HOME", str(assets_home))
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)

    output_dir = tmp_path / "db"
    exit_code = cli.main(
        [
            "build-db",
            "--offline",
            "--output-dir",
            str(output_dir),
        ]
    )

    assert exit_code == 0
    assert (output_dir / "marker_db.parquet").exists()
    assert (output_dir / "marker_db.sqlite").exists()
