"""LLM-powered annotation engine for GPT Cell Annotator."""

from __future__ import annotations

import json
import logging
import time
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any, cast

from jsonschema import Draft202012Validator

try:  # structlog is optional in certain environments (e.g., tests)
    from structlog.contextvars import get_contextvars
except ImportError:  # pragma: no cover - fallback when structlog not available

    def get_contextvars() -> dict[str, Any]:
        return {}


from openai import OpenAI

from backend.llm import prompts
from backend.llm.retrieval import retrieve_candidates
from backend.util.gene_normalization import get_gene_normalizer
from config.settings import Settings, get_settings

logger = logging.getLogger("gpt_cell_annotator.llm")


class AnnotationError(RuntimeError):
    """Raised when the LLM returns an invalid response."""


class SchemaValidationError(AnnotationError):
    """Raised when LLM output fails JSON Schema validation."""


DEFAULT_MARKER_KB: dict[str, dict[str, Any]] = {
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

    def __init__(self, knowledge_base: dict[str, dict[str, Any]] | None = None) -> None:
        self.knowledge_base = knowledge_base or DEFAULT_MARKER_KB

    def annotate_cluster(
        self,
        cluster_payload: dict[str, Any],
        dataset_context: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
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
            alternatives: list[dict[str, str]] = []
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
        clusters: Iterable[dict[str, Any]],
        dataset_context: dict[str, Any] | None = None,
    ) -> dict[str, dict[str, Any]]:
        results: dict[str, dict[str, Any]] = {}
        for cluster in clusters:
            cluster_id = str(cluster.get("cluster_id", "unknown"))
            results[cluster_id] = self.annotate_cluster(cluster, dataset_context)
        return results

    @staticmethod
    def _normalise_markers(markers: Sequence[str]) -> set[str]:
        return {
            marker.upper().strip()
            for marker in markers
            if isinstance(marker, str) and marker.strip()
        }

    def _rank_labels(
        self,
        markers: set[str],
        dataset_context: dict[str, Any] | None,
    ) -> list[tuple[int, str, set[str]]]:
        ranked: list[tuple[int, str, set[str]]] = []
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

    def _build_alternatives(
        self,
        ranked: Iterable[tuple[int, str, set[str]]],
    ) -> list[dict[str, str]]:
        alternatives: list[dict[str, str]] = []
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
        settings: Settings | None = None,
        *,
        client: OpenAI | None = None,
        mock_backend: MockAnnotator | None = None,
        force_mock: bool = False,
    ) -> None:
        self.settings = settings or get_settings()
        self._mock_backend = mock_backend or MockAnnotator()
        self._client: OpenAI | None
        self._gene_normalizer = get_gene_normalizer()
        self._primary_species = self._gene_normalizer.primary_species
        self._force_mock = force_mock

        if force_mock:
            self._client = None
        elif client is not None:
            self._client = client
        elif self.settings.openai_api_key:
            self._client = OpenAI(api_key=self.settings.openai_api_key)
        else:
            self._client = None

        self._mode = "live" if self._client is not None else "mock"
        if self._mode == "mock":
            if self._force_mock:
                logger.info(
                    "Running Annotator in offline mock mode; external LLM calls are disabled."
                )
            else:
                logger.warning(
                    "OPENAI_API_KEY not configured; using heuristic MockAnnotator. "
                    "Outputs are suitable for demos only."
                )

        self._min_interval = (
            0.0
            if self.settings.openai_requests_per_minute <= 0
            else 60.0 / self.settings.openai_requests_per_minute
        )
        self._last_call_ts: float | None = None
        self._annotation_schema: dict[str, Any]
        self._annotation_validator: Draft202012Validator
        self._batch_schema: dict[str, Any]
        self._batch_validator: Draft202012Validator
        self._load_schemas()

    @property
    def llm_mode(self) -> str:
        return self._mode

    # Public API -----------------------------------------------------------------

    def _load_schemas(self) -> None:
        schema_path = (
            Path(__file__).resolve().parents[2] / "schemas" / "annotation_result.schema.json"
        )
        try:
            with schema_path.open("r", encoding="utf-8") as fh:
                dataset_schema = json.load(fh)
        except FileNotFoundError:
            logger.error(
                "Annotation schema not found at %s; using permissive validator.",
                schema_path,
            )
            dataset_schema = {}

        annotation_schema = (
            dataset_schema.get("properties", {})
            .get("clusters", {})
            .get("items", {})
            .get("properties", {})
            .get("annotation")
        )
        if not annotation_schema:
            logger.warning(
                "Annotation schema missing from %s; validation will allow any object.", schema_path
            )
            annotation_schema = {"type": "object"}

        self._annotation_schema = annotation_schema
        self._annotation_validator = Draft202012Validator(self._annotation_schema)
        self._batch_schema = {
            "type": "object",
            "additionalProperties": self._annotation_schema,
        }
        self._batch_validator = Draft202012Validator(self._batch_schema)

    def annotate_cluster(
        self,
        cluster_payload: dict[str, Any],
        dataset_context: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Annotate a single cluster and return the parsed JSON response."""

        enriched_payload = self._prepare_cluster_payload(
            cluster_payload,
            dataset_context,
        )

        if self._mode == "mock" or self._client is None:
            return self._mock_backend.annotate_cluster(
                self._mock_payload(enriched_payload),
                dataset_context,
            )

        self._log_request("annotate_cluster", enriched_payload, dataset_context)
        messages = [
            {"role": "system", "content": prompts.build_system_prompt()},
            {
                "role": "user",
                "content": prompts.build_single_cluster_prompt(
                    enriched_payload,
                    dataset_context,
                ),
            },
        ]
        try:
            raw = self._call_llm(
                messages,
                schema=self._annotation_schema,
                schema_name="gpt_cell_annotator_annotation",
            )
            result = self._parse_json(raw, context="annotate_cluster")
            self._validate_annotation_payload(result)
            self._log_response("annotate_cluster", result)
            mapping_notes = enriched_payload.get("mapping_notes", [])
            if mapping_notes:
                metadata = result.setdefault("metadata", {})
                metadata["mapping_notes"] = mapping_notes
                metadata["canonical_markers"] = enriched_payload.get("canonical_markers")
                metadata["target_species"] = enriched_payload.get("_target_species")
                metadata["mapping_applied"] = enriched_payload.get("_mapping_applied")
                metadata["original_markers"] = enriched_payload.get("original_markers")
            return result
        except (AnnotationError, SchemaValidationError) as exc:
            logger.warning("annotate_cluster falling back to heuristics: %s", exc)
            fallback = self._mock_backend.annotate_cluster(
                self._mock_payload(enriched_payload),
                dataset_context,
            )
            self._append_warning(fallback, f"Mock annotator used due to LLM error: {exc}")
            return fallback

    def annotate_batch(
        self,
        clusters: Iterable[dict[str, Any]],
        dataset_context: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Annotate multiple clusters in one call."""

        enriched_clusters = [
            self._prepare_cluster_payload(cluster, dataset_context) for cluster in clusters
        ]
        if self._mode == "mock" or self._client is None:
            return self._mock_backend.annotate_batch(
                [self._mock_payload(cluster) for cluster in enriched_clusters],
                dataset_context,
            )
        self._log_request("annotate_batch", enriched_clusters, dataset_context)
        messages = [
            {"role": "system", "content": prompts.build_system_prompt()},
            {
                "role": "user",
                "content": prompts.build_batch_prompt(enriched_clusters, dataset_context),
            },
        ]
        try:
            raw = self._call_llm(
                messages,
                schema=self._batch_schema,
                schema_name="gpt_cell_annotator_batch",
            )
            result = self._parse_json(raw, context="annotate_batch")
            self._validate_batch_payload(result)
            self._log_response("annotate_batch", result)
            for original_cluster, enriched in zip(clusters, enriched_clusters, strict=True):
                cluster_id = str(original_cluster.get("cluster_id"))
                metadata = result.setdefault(cluster_id, {}).setdefault("metadata", {})
                if enriched.get("mapping_notes"):
                    metadata["mapping_notes"] = enriched.get("mapping_notes")
                if enriched.get("canonical_markers"):
                    metadata["canonical_markers"] = enriched.get("canonical_markers")
                metadata["target_species"] = enriched.get("_target_species")
                metadata["mapping_applied"] = enriched.get("_mapping_applied")
                metadata["original_markers"] = enriched.get("original_markers")
            return result
        except (AnnotationError, SchemaValidationError) as exc:
            logger.warning("annotate_batch falling back to heuristics: %s", exc)
            fallback = self._mock_backend.annotate_batch(
                [self._mock_payload(cluster) for cluster in enriched_clusters],
                dataset_context,
            )
            for cluster_result in fallback.values():
                self._append_warning(cluster_result, f"Mock annotator used due to LLM error: {exc}")
            return fallback

    # Internal helpers -----------------------------------------------------------

    def _parse_json(self, raw: str, *, context: str) -> dict[str, Any]:
        try:
            data = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise AnnotationError(f"{context} returned non-JSON payload: {raw!r}") from exc
        if not isinstance(data, dict):
            raise AnnotationError(f"{context} expected JSON object, received {type(data)}")
        return data

    def _validate_annotation_payload(self, payload: dict[str, Any]) -> None:
        if not isinstance(payload, dict):
            raise SchemaValidationError("Annotation payload must be an object.")
        errors = sorted(self._annotation_validator.iter_errors(payload), key=lambda err: err.path)
        if errors:
            first = errors[0]
            context = get_contextvars()
            trace_id = context.get("trace_id")
            logger.error(
                "annotation.schema_validation_failed",
                extra={
                    "error": first.message,
                    "path": list(first.path),
                    "schema_path": list(first.schema_path),
                    "trace_id": trace_id,
                },
            )
            raise SchemaValidationError(first.message)

    def _validate_batch_payload(self, payload: dict[str, Any]) -> None:
        if not isinstance(payload, dict):
            raise SchemaValidationError("Batch payload must be an object mapping cluster IDs.")
        errors = sorted(self._batch_validator.iter_errors(payload), key=lambda err: err.path)
        if errors:
            first = errors[0]
            context = get_contextvars()
            trace_id = context.get("trace_id")
            logger.error(
                "annotation.batch_schema_validation_failed",
                extra={
                    "error": first.message,
                    "path": list(first.path),
                    "schema_path": list(first.schema_path),
                    "trace_id": trace_id,
                },
            )
            raise SchemaValidationError(first.message)

    @staticmethod
    def _append_warning(annotation: dict[str, Any], message: str) -> None:
        warnings = annotation.get("warnings")
        if isinstance(warnings, list):
            warnings.append(message)
        else:
            annotation["warnings"] = [message]

    @staticmethod
    def _mock_payload(cluster_payload: dict[str, Any]) -> dict[str, Any]:
        mock_payload = dict(cluster_payload)
        canonical = cluster_payload.get("canonical_markers")
        if canonical:
            mock_payload["markers"] = canonical
        return mock_payload

    def _prepare_cluster_payload(
        self,
        cluster_payload: dict[str, Any],
        dataset_context: dict[str, Any] | None,
    ) -> dict[str, Any]:
        payload = dict(cluster_payload)
        context = dataset_context or {}
        species = context.get("species") if isinstance(context, dict) else None
        tissue = context.get("tissue") if isinstance(context, dict) else None

        original_markers = [
            marker
            for marker in (payload.get("markers") or [])
            if isinstance(marker, str) and marker.strip()
        ]
        payload.setdefault("markers", original_markers)
        payload["original_markers"] = list(original_markers)

        canonical_markers, mapping_notes = self._gene_normalizer.map_to_primary(
            original_markers,
            species,
        )
        mapping_applied = bool(mapping_notes) or (
            species and species.lower() != self._primary_species.lower()
        )
        target_species = (
            self._primary_species
            if species and species.lower() != self._primary_species.lower()
            else (species or self._primary_species)
        )

        payload["canonical_markers"] = canonical_markers or original_markers
        payload["mapping_notes"] = mapping_notes
        payload["_mapping_notes"] = mapping_notes
        payload["_mapping_applied"] = mapping_applied
        payload["_target_species"] = target_species

        markers_for_retrieval = payload["canonical_markers"] or original_markers

        if self.settings.rag_enabled:
            candidates = retrieve_candidates(
                markers_for_retrieval,
                species=target_species,
                tissue=tissue,
            )
            if candidates:
                payload["retrieved_candidates"] = [
                    {
                        "cell_type": candidate.cell_type,
                        "ontology_id": candidate.ontology_id,
                        "overlap": candidate.overlap,
                        "supporting_markers": candidate.supporting_markers,
                        "tissue_counts": candidate.tissue_counts,
                    }
                    for candidate in candidates
                ]
        return payload

    def _log_request(
        self,
        operation: str,
        payload: Iterable[dict[str, Any]] | dict[str, Any],
        dataset_context: dict[str, Any] | None,
    ) -> None:
        if isinstance(payload, dict):
            marker_counts = len(payload.get("markers") or [])
            cluster_count = 1
        else:
            payload_list = list(payload)
            marker_counts = sum(len(item.get("markers") or []) for item in payload_list)
            cluster_count = len(payload_list)
        context_keys = sorted((dataset_context or {}).keys())
        logger.info(
            "llm.request",
            extra={
                "operation": operation,
                "cluster_count": cluster_count,
                "marker_count": marker_counts,
                "context_keys": context_keys,
            },
        )

    def _log_response(self, operation: str, payload: dict[str, Any]) -> None:
        logger.info(
            "llm.response",
            extra={
                "operation": operation,
                "keys": sorted(payload.keys()),
            },
        )

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

    def _call_llm(
        self,
        messages: list[dict[str, str]],
        *,
        schema: dict[str, Any],
        schema_name: str,
    ) -> str:
        if self._client is None:
            raise AnnotationError("LLM client unavailable; mock mode should be used instead.")

        retries = max(self.settings.openai_retry_attempts, 1)
        for attempt in range(1, retries + 1):
            self._enforce_rate_limit()
            try:
                response_format: Any = {
                    "type": "json_schema",
                    "json_schema": {
                        "name": schema_name,
                        "schema": schema,
                    },
                }
                response = self._client.chat.completions.create(
                    model=self.settings.openai_model,
                    messages=cast(Any, messages),
                    temperature=self.settings.openai_temperature,
                    max_tokens=self.settings.openai_max_tokens,
                    response_format=cast(Any, response_format),
                )
                content = response.choices[0].message.content
                if content is None:
                    raise AnnotationError("LLM returned empty content")
                if isinstance(content, list):
                    raise AnnotationError("LLM returned non-text content")
                return str(content).strip()
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


__all__ = ["AnnotationError", "Annotator", "MockAnnotator", "SchemaValidationError"]
