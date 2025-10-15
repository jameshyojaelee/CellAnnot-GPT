from __future__ import annotations

import json
from types import SimpleNamespace
from typing import Any

import pytest

from backend.llm.annotator import Annotator
from config.settings import Settings


class FakeCompletions:
    def __init__(self, responses: list[Any]) -> None:
        self._responses = responses
        self.calls: list[dict[str, Any]] = []

    def create(self, **kwargs: Any) -> Any:
        self.calls.append(kwargs)
        response = self._responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


class FakeClient:
    def __init__(self, completions: FakeCompletions) -> None:
        self.chat = SimpleNamespace(completions=completions)


class SleepRecorder:
    def __init__(self) -> None:
        self.calls: list[float] = []

    def __call__(self, seconds: float) -> None:
        self.calls.append(seconds)


def build_choice(content: str) -> Any:
    return SimpleNamespace(choices=[SimpleNamespace(message=SimpleNamespace(content=content))])


def test_annotate_cluster_parses_json(monkeypatch: pytest.MonkeyPatch) -> None:
    expected = {
        "primary_label": "B cell",
        "confidence": "High",
        "rationale": "Markers support B cell",
    }
    fake_responses = [build_choice(json.dumps(expected))]
    completions = FakeCompletions(fake_responses)
    client = FakeClient(completions)
    settings = Settings(
        openai_api_key="test-key",
        openai_requests_per_minute=0,
        openai_retry_attempts=1,
    )
    annotator = Annotator(settings=settings, client=client)  # type: ignore[arg-type]
    assert annotator.llm_mode == "live"
    recorder = SleepRecorder()
    monkeypatch.setattr(annotator, "_sleep", recorder)  # type: ignore[arg-type]

    payload = {"cluster_id": "0", "markers": ["MS4A1"]}
    result = annotator.annotate_cluster(payload)

    assert result == expected
    assert len(completions.calls) == 1
    assert recorder.calls == []
    assert completions.calls[0]["response_format"]["type"] == "json_schema"


def test_annotate_batch_retries_on_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    first_error = RuntimeError("rate limit")
    expected = {
        "0": {
            "primary_label": "T cell",
            "confidence": "Medium",
            "rationale": "Markers align with canonical T cell profile",
        }
    }
    fake_responses = [first_error, build_choice(json.dumps(expected))]
    completions = FakeCompletions(fake_responses)
    client = FakeClient(completions)
    settings = Settings(
        openai_api_key="test-key",
        openai_requests_per_minute=0,
        openai_retry_attempts=2,
        openai_retry_backoff_seconds=0.1,
    )
    annotator = Annotator(settings=settings, client=client)  # type: ignore[arg-type]
    assert annotator.llm_mode == "live"
    recorder = SleepRecorder()
    monkeypatch.setattr(annotator, "_sleep", recorder)  # type: ignore[arg-type]

    result = annotator.annotate_batch([{"cluster_id": "0", "markers": ["CD3E"]}])

    assert result == expected
    assert len(completions.calls) == 2  # retried once
    assert pytest.approx(recorder.calls[0], rel=1e-6) == 0.1
    assert completions.calls[1]["response_format"]["type"] == "json_schema"


def test_invalid_json_triggers_mock_fallback() -> None:
    completions = FakeCompletions([build_choice("not json")])
    client = FakeClient(completions)
    settings = Settings(openai_api_key="test-key", openai_requests_per_minute=0)
    annotator = Annotator(settings=settings, client=client)  # type: ignore[arg-type]
    assert annotator.llm_mode == "live"

    result = annotator.annotate_cluster({"cluster_id": "1", "markers": []})

    assert result["primary_label"] in {"Unknown or Novel", "B cell", "T cell"}
    assert "warnings" in result
    assert any("Mock annotator used" in warning for warning in result["warnings"])


def test_mock_mode_when_api_key_missing():
    settings = Settings(openai_api_key="", openai_requests_per_minute=0)
    annotator = Annotator(settings=settings)

    assert annotator.llm_mode == "mock"
    result = annotator.annotate_cluster({"cluster_id": "0", "markers": ["MS4A1", "CD79A"]})
    assert result["primary_label"] == "B cell"
    assert "rationale" in result


def test_schema_validation_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    bad_payload = {"confidence": "High", "rationale": "Missing primary label"}
    completions = FakeCompletions([build_choice(json.dumps(bad_payload))])
    client = FakeClient(completions)
    settings = Settings(openai_api_key="test-key", openai_requests_per_minute=0)
    annotator = Annotator(settings=settings, client=client)  # type: ignore[arg-type]

    result = annotator.annotate_cluster({"cluster_id": "7", "markers": ["MS4A1"]})

    assert result["primary_label"] == "B cell"
    assert "warnings" in result
    assert any("LLM error" in warning for warning in result["warnings"])


def test_batch_schema_validation_fallback() -> None:
    bad_batch = {"0": {"rationale": "Missing confidence and primary_label"}}
    completions = FakeCompletions([build_choice(json.dumps(bad_batch))])
    client = FakeClient(completions)
    settings = Settings(openai_api_key="test-key", openai_requests_per_minute=0)
    annotator = Annotator(settings=settings, client=client)  # type: ignore[arg-type]

    result = annotator.annotate_batch([{"cluster_id": "0", "markers": ["CD3E"]}])

    assert isinstance(result["0"]["primary_label"], str)
    assert "warnings" in result["0"]
    assert any("LLM error" in warning for warning in result["0"]["warnings"])
