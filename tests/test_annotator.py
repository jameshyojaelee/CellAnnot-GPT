from __future__ import annotations

import json
from types import SimpleNamespace
from typing import Any, Dict, List

import pytest

from backend.llm.annotator import Annotator, AnnotationError
from config.settings import Settings


class FakeCompletions:
    def __init__(self, responses: List[Any]) -> None:
        self._responses = responses
        self.calls: List[Dict[str, Any]] = []

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
        self.calls: List[float] = []

    def __call__(self, seconds: float) -> None:
        self.calls.append(seconds)


def build_choice(content: str) -> Any:
    return SimpleNamespace(choices=[SimpleNamespace(message=SimpleNamespace(content=content))])


def test_annotate_cluster_parses_json(monkeypatch: pytest.MonkeyPatch) -> None:
    expected = {"primary_label": "B cell"}
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

    result = annotator.annotate_cluster({"cluster_id": "0", "markers": ["MS4A1"]})

    assert result == expected
    assert len(completions.calls) == 1
    assert recorder.calls == []


def test_annotate_batch_retries_on_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    first_error = RuntimeError("rate limit")
    expected = {"Cluster 0": {"primary_label": "T cell"}}
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


def test_invalid_json_raises_annotation_error() -> None:
    completions = FakeCompletions([build_choice("not json")])
    client = FakeClient(completions)
    settings = Settings(openai_api_key="test-key", openai_requests_per_minute=0)
    annotator = Annotator(settings=settings, client=client)  # type: ignore[arg-type]
    assert annotator.llm_mode == "live"

    with pytest.raises(AnnotationError):
        annotator.annotate_cluster({"cluster_id": "1", "markers": []})


def test_mock_mode_when_api_key_missing():
    settings = Settings(openai_api_key="", openai_requests_per_minute=0)
    annotator = Annotator(settings=settings)

    assert annotator.llm_mode == "mock"
    result = annotator.annotate_cluster({"cluster_id": "0", "markers": ["MS4A1", "CD79A"]})
    assert result["primary_label"] == "B cell"
    assert "rationale" in result
