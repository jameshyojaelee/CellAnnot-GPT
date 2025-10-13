"""LLM-powered annotation engine for CellAnnot-GPT."""

from __future__ import annotations

import json
import time
from typing import Any, Dict, Iterable, List, Optional

from openai import OpenAI

from backend.llm import prompts
from config.settings import Settings, get_settings


class AnnotationError(RuntimeError):
    """Raised when the LLM returns an invalid response."""


class Annotator:
    """High-level orchestrator that wraps prompt construction and OpenAI calls."""

    def __init__(
        self,
        settings: Optional[Settings] = None,
        *,
        client: Optional[OpenAI] = None,
    ) -> None:
        self.settings = settings or get_settings()
        self._client = client or OpenAI(api_key=self.settings.openai_api_key or None)
        self._min_interval = (
            0.0
            if self.settings.openai_requests_per_minute <= 0
            else 60.0 / self.settings.openai_requests_per_minute
        )
        self._last_call_ts: Optional[float] = None

    # Public API -----------------------------------------------------------------

    def annotate_cluster(
        self,
        cluster_payload: Dict[str, Any],
        dataset_context: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Annotate a single cluster and return the parsed JSON response."""

        messages = [
            {"role": "system", "content": prompts.build_system_prompt()},
            {
                "role": "user",
                "content": prompts.build_single_cluster_prompt(cluster_payload, dataset_context),
            },
        ]
        raw = self._call_llm(messages)
        return self._parse_json(raw, context="annotate_cluster")

    def annotate_batch(
        self,
        clusters: Iterable[Dict[str, Any]],
        dataset_context: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Annotate multiple clusters in one call."""

        messages = [
            {"role": "system", "content": prompts.build_system_prompt()},
            {
                "role": "user",
                "content": prompts.build_batch_prompt(clusters, dataset_context),
            },
        ]
        raw = self._call_llm(messages)
        return self._parse_json(raw, context="annotate_batch")

    # Internal helpers -----------------------------------------------------------

    def _parse_json(self, raw: str, *, context: str) -> Dict[str, Any]:
        try:
            data = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise AnnotationError(f"{context} returned non-JSON payload: {raw!r}") from exc
        if not isinstance(data, dict):
            raise AnnotationError(f"{context} expected JSON object, received {type(data)}")
        return data

    def _enforce_rate_limit(self) -> None:
        if self._min_interval <= 0:
            return
        now = time.monotonic()
        if self._last_call_ts is None:
            self._last_call_ts = now
            return
        elapsed = now - self._last_call_ts
        wait_for = self._min_interval - elapsed
        if wait_for > 0:
            self._sleep(wait_for)
        self._last_call_ts = time.monotonic()

    def _call_llm(self, messages: List[Dict[str, str]]) -> str:
        retries = max(self.settings.openai_retry_attempts, 1)
        for attempt in range(1, retries + 1):
            self._enforce_rate_limit()
            try:
                response = self._client.chat.completions.create(
                    model=self.settings.openai_model,
                    messages=messages,
                    temperature=self.settings.openai_temperature,
                    max_tokens=self.settings.openai_max_tokens,
                )
                content = response.choices[0].message.content  # type: ignore[attr-defined]
                if not content:
                    raise AnnotationError("LLM returned empty content")
                return content.strip()
            except AnnotationError:
                raise
            except Exception as exc:  # pragma: no cover - fallback path tested via mocks
                if attempt == retries:
                    raise
                backoff = self.settings.openai_retry_backoff_seconds * attempt
                self._sleep(backoff)
                continue

        raise AnnotationError("Failed to obtain response from LLM")  # pragma: no cover

    def _sleep(self, seconds: float) -> None:
        time.sleep(seconds)


__all__ = ["Annotator", "AnnotationError"]
