"""LLM-powered annotation engine for CellAnnot-GPT."""

from __future__ import annotations

import json
import logging
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from openai import OpenAI

from backend.llm import prompts
from config.settings import Settings, get_settings

logger = logging.getLogger("cellannot.llm")


class AnnotationError(RuntimeError):
    """Raised when the LLM returns an invalid response."""


DEFAULT_MARKER_KB: Dict[str, Dict[str, Any]] = {
    "B cell": {
        "markers": {"MS4A1", "CD79A", "CD74", "CD19"},
        "ontology_id": "CL:0000236",
    },
    "CD4 T cell": {
        "markers": {"CD3D", "CD3E", "CD4", "IL7R"},
        "ontology_id": "CL:0000624",
    },
    "CD8 T cell": {
        "markers": {"CD3D", "CD3E", "CD8A", "GZMB"},
        "ontology_id": "CL:0000625",
    },
    "NK cell": {
        "markers": {"NKG7", "GNLY", "PRF1", "KLRD1"},
        "ontology_id": "CL:0000623",
    },
    "Monocyte": {
        "markers": {"LYZ", "S100A9", "LGALS3", "FCGR3A", "MS4A7"},
        "ontology_id": "CL:0000576",
    },
    "Platelet": {
        "markers": {"PPBP", "PF4", "NRGN"},
        "ontology_id": "CL:0000233",
    },
    "Erythrocyte": {
        "markers": {"HBB", "HBA1", "AHSP"},
        "ontology_id": "CL:0000232",
    },
}


class MockAnnotator:
    """Heuristic annotator used when external LLM access is unavailable."""

    def __init__(self, knowledge_base: Optional[Dict[str, Dict[str, Any]]] = None) -> None:
        self.knowledge_base = knowledge_base or DEFAULT_MARKER_KB

    def annotate_cluster(
        self,
        cluster_payload: Dict[str, Any],
        dataset_context: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        markers = self._normalise_markers(cluster_payload.get("markers", []))
        ranked = self._rank_labels(markers, dataset_context)

        if ranked:
            best_score, best_label, shared = ranked[0]
        else:
            best_score, best_label, shared = (0, "Unknown or Novel", set())

        if best_score == 0:
            primary_label = "Unknown or Novel"
            ontology_id = None
            rationale = "No strong marker overlap; reporting as Unknown/Novel."
            confidence = "Low"
            caveats = "Heuristic fallback due to missing LLM output."
            alternatives: List[Dict[str, str]] = []
        else:
            primary_label = best_label
            ontology_id = self.knowledge_base[best_label].get("ontology_id")
            rationale = f"Matched markers: {', '.join(sorted(shared))}"
            confidence = self._confidence_from_score(best_score)
            caveats = None
            alternatives = self._build_alternatives(ranked[1:])

        return {
            "primary_label": primary_label,
            "ontology_id": ontology_id,
            "confidence": confidence,
            "rationale": rationale,
            "alternatives": alternatives,
            "caveats": caveats,
        }

    def annotate_batch(
        self,
        clusters: Iterable[Dict[str, Any]],
        dataset_context: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Dict[str, Any]]:
        results: Dict[str, Dict[str, Any]] = {}
        for cluster in clusters:
            cluster_id = str(cluster.get("cluster_id", "unknown"))
            results[cluster_id] = self.annotate_cluster(cluster, dataset_context)
        return results

    @staticmethod
    def _normalise_markers(markers: Sequence[str]) -> Set[str]:
        return {marker.upper().strip() for marker in markers if isinstance(marker, str) and marker.strip()}

    def _rank_labels(
        self,
        markers: Set[str],
        dataset_context: Optional[Dict[str, Any]],
    ) -> List[Tuple[int, str, Set[str]]]:
        ranked: List[Tuple[int, str, Set[str]]] = []
        for label, info in self.knowledge_base.items():
            kb_markers = info.get("markers", set())
            shared = markers & kb_markers
            if shared:
                ranked.append((len(shared), label, shared))
        ranked.sort(key=lambda item: (item[0], item[1]), reverse=True)
        return ranked

    @staticmethod
    def _confidence_from_score(score: int) -> str:
        if score >= 3:
            return "High"
        if score == 2:
            return "Medium"
        return "Low"

    def _build_alternatives(self, ranked: Iterable[Tuple[int, str, Set[str]]]) -> List[Dict[str, str]]:
        alternatives: List[Dict[str, str]] = []
        for score, label, shared in ranked:
            if score == 0 or len(alternatives) >= 2:
                break
            alternatives.append(
                {
                    "label": label,
                    "reason": f"Shared markers: {', '.join(sorted(shared))}",
                }
            )
        return alternatives


class Annotator:
    """High-level orchestrator that wraps prompt construction and OpenAI calls."""

    def __init__(
        self,
        settings: Optional[Settings] = None,
        *,
        client: Optional[OpenAI] = None,
        mock_backend: Optional[MockAnnotator] = None,
    ) -> None:
        self.settings = settings or get_settings()
        self._mock_backend = mock_backend or MockAnnotator()
        self._client: Optional[OpenAI]

        if client is not None:
            self._client = client
        elif self.settings.openai_api_key:
            self._client = OpenAI(api_key=self.settings.openai_api_key)
        else:
            self._client = None

        self._mode = "live" if self._client is not None else "mock"
        if self._mode == "mock":
            logger.warning(
                "OPENAI_API_KEY not configured; using heuristic MockAnnotator. "
                "Outputs are suitable for demos only."
            )

        self._min_interval = (
            0.0
            if self.settings.openai_requests_per_minute <= 0
            else 60.0 / self.settings.openai_requests_per_minute
        )
        self._last_call_ts: Optional[float] = None

    @property
    def llm_mode(self) -> str:
        return self._mode

    # Public API -----------------------------------------------------------------

    def annotate_cluster(
        self,
        cluster_payload: Dict[str, Any],
        dataset_context: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Annotate a single cluster and return the parsed JSON response."""

        if self._mode == "mock" or self._client is None:
            return self._mock_backend.annotate_cluster(cluster_payload, dataset_context)

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

        if self._mode == "mock" or self._client is None:
            return self._mock_backend.annotate_batch(clusters, dataset_context)

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
        if self._client is None:
            raise AnnotationError("LLM client unavailable; mock mode should be used instead.")

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
                logger.warning("LLM call failed (attempt %s/%s): %s", attempt, retries, exc)
                if attempt == retries:
                    raise
                backoff = self.settings.openai_retry_backoff_seconds * attempt
                self._sleep(backoff)
                continue

        raise AnnotationError("Failed to obtain response from LLM")  # pragma: no cover

    def _sleep(self, seconds: float) -> None:
        time.sleep(seconds)


__all__ = ["Annotator", "AnnotationError", "MockAnnotator"]
