"""Scanpy integration helpers, annotation APIs, and CLI entrypoints."""

from __future__ import annotations

import argparse
import asyncio
import json
import logging
import os
import sys
import time
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import (
    Any,
    Iterable,
    Mapping,
    NamedTuple,
    Protocol,
    Sequence,
    TYPE_CHECKING,
)
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData

from backend.llm.annotator import Annotator
from backend.validation.crosscheck import crosscheck_batch
from backend.validation.report import DatasetReport, build_structured_report
from config.settings import Settings, get_settings
from gpt_cell_annotator import assets

try:  # Optional dependency used for ranking and Loom IO.
    import scanpy as sc
except ImportError:  # pragma: no cover - handled via guards in runtime paths.
    sc = None

try:  # Prefer structlog for telemetry; fall back to stdlib logging otherwise.
    import structlog

    _logger = structlog.get_logger("gpt_cell_annotator.scanpy")

    def _log(event: str, /, **fields: Any) -> None:
        _logger.info(event, **fields)


except ImportError:  # pragma: no cover - structlog is an optional extra at runtime.
    structlog = None  # type: ignore
    _logger = logging.getLogger("gpt_cell_annotator.scanpy")

    def _log(event: str, /, **fields: Any) -> None:
        _logger.info("%s %s", event, fields)


if TYPE_CHECKING:  # pragma: no cover - typing only import
    from backend.cache.redis_cache import RedisAnnotationCache

try:  # Optional metrics sink for extras[gca,api]
    from prometheus_client import Counter, Histogram
except ImportError:  # pragma: no cover - optional dependency
    Counter = Histogram = None  # type: ignore


MARKER_DB_COLUMNS = [
    "source",
    "cell_type",
    "ontology_id",
    "gene_symbol",
    "species",
    "tissue",
    "evidence",
    "reference",
    "evidence_score",
]

CACHE_VERSION = "scanpy-v2"
MARKER_DB_ENV_VAR = "GCA_MARKER_DB_PATH"
CACHE_DIR_ENV_VAR = "GCA_CACHE_DIR"
REQUEST_ID_ENV_VAR = "GCA_REQUEST_ID"

SPECIES_PRESETS: dict[str, dict[str, str | None]] = {
    "human_pbmc": {"species": "Homo sapiens", "tissue": "Peripheral blood"},
    "mouse_pbmc": {"species": "Mus musculus", "tissue": "Peripheral blood"},
    "mouse_brain": {"species": "Mus musculus", "tissue": "Brain"},
}

if Counter is not None:
    _BATCH_COUNTER = Counter(
        "gca_scanpy_batches_total",
        "Number of annotation batches executed through the Scanpy integration.",
        ("mode",),
    )
    _CLUSTER_COUNTER = Counter(
        "gca_scanpy_clusters_total",
        "Number of clusters annotated through the Scanpy integration.",
        ("mode",),
    )
    _BATCH_DURATION = Histogram(
        "gca_scanpy_batch_duration_seconds",
        "Duration of annotation batches in seconds.",
        ("mode",),
    )
else:  # pragma: no cover - metrics optional
    _BATCH_COUNTER = None
    _CLUSTER_COUNTER = None
    _BATCH_DURATION = None


class AnnotationCacheProtocol(Protocol):
    """Minimal async cache interface used by the batching workflow."""

    async def get(self, payload: dict[str, Any]) -> dict[str, Any] | None:  # pragma: no cover - protocol
        ...

    async def set(self, payload: dict[str, Any], value: dict[str, Any]) -> None:  # pragma: no cover - protocol
        ...


class DiskAnnotationCache:
    """Disk-backed cache for Scanpy annotations (offline-friendly)."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)

    def _path_for(self, payload: dict[str, Any]) -> Path:
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        digest = sha256(canonical.encode("utf-8")).hexdigest()
        return self.root / f"{digest}.json"

    async def get(self, payload: dict[str, Any]) -> dict[str, Any] | None:
        path = self._path_for(payload)
        if not path.exists():
            return None
        return await asyncio.to_thread(self._read, path)

    async def set(self, payload: dict[str, Any], value: dict[str, Any]) -> None:
        path = self._path_for(payload)
        await asyncio.to_thread(self._write, path, value)

    @staticmethod
    def _read(path: Path) -> dict[str, Any] | None:
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None

    @staticmethod
    def _write(path: Path, value: dict[str, Any]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value), encoding="utf-8")


@dataclass(slots=True)
class GuardrailConfig:
    """Overrides for guardrail behaviour and confidence calibration."""

    min_marker_overlap: int | None = None
    force_unknown_on_fail: bool | None = None
    confidence_overlap_medium: int | None = None
    confidence_overlap_high: int | None = None

    def as_update(self) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if self.min_marker_overlap is not None:
            payload["validation_min_marker_overlap"] = max(0, self.min_marker_overlap)
        if self.force_unknown_on_fail is not None:
            payload["validation_force_unknown_on_fail"] = bool(self.force_unknown_on_fail)
        if self.confidence_overlap_medium is not None:
            payload["confidence_overlap_medium"] = max(0, self.confidence_overlap_medium)
        if self.confidence_overlap_high is not None:
            payload["confidence_overlap_high"] = max(
                self.confidence_overlap_medium or 0,
                self.confidence_overlap_high,
            )
        return payload


@dataclass(slots=True)
class BatchOptions:
    """Controls batching characteristics for annotation requests."""

    size: int = 32
    concurrency: int = 1

    def normalized(self, total: int) -> BatchOptions:
        size = max(1, self.size or total)
        concurrency = max(1, self.concurrency)
        if total < size:
            size = total
        if total < concurrency:
            concurrency = total
        return BatchOptions(size=size, concurrency=concurrency)


class ScanpyAnnotationResult(NamedTuple):
    """Structured response returned from annotate_anndata."""

    adata: AnnData
    report: DatasetReport
    annotations: list[dict[str, Any]]
    stats: dict[str, Any]

    @property
    def report_dict(self) -> dict[str, Any]:
        return self.report.model_dump()


class MarkerAnnotationResult(NamedTuple):
    """Result returned when annotating pre-computed marker collections."""

    annotations: list[dict[str, Any]]
    report: DatasetReport
    stats: dict[str, Any]

    @property
    def report_dict(self) -> dict[str, Any]:
        return self.report.model_dump()


def _default_request_id(user_supplied: str | None) -> str:
    env_override = os.environ.get(REQUEST_ID_ENV_VAR)
    if user_supplied:
        return user_supplied
    if env_override:
        return env_override
    return str(uuid4())


async def _run_in_thread(func: Any, *args: Any, **kwargs: Any) -> Any:
    return await asyncio.to_thread(lambda: func(*args, **kwargs))


def _ensure_rankings(
    adata: AnnData,
    cluster_key: str,
    *,
    top_n_markers: int,
    method: str = "wilcoxon",
) -> None:
    """Compute rank_genes_groups if missing or incompatible with cluster key."""

    rankings = adata.uns.get("rank_genes_groups")
    params = (rankings or {}).get("params", {})
    if rankings is not None and params.get("groupby") == cluster_key:
        return

    if sc is None:
        raise ImportError(
            "scanpy is required to compute marker rankings. Install the "
            'extra with `pip install "gpt-cell-annotator[scanpy]"` or '
            "precompute `rank_genes_groups` before calling annotate_anndata."
        )

    sc.tl.rank_genes_groups(adata, groupby=cluster_key, n_genes=top_n_markers, method=method)


def _unique_ordered(markers: Iterable[Any]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for marker in markers:
        if not isinstance(marker, str):
            continue
        normalised = marker.strip()
        if not normalised:
            continue
        upper = normalised.upper()
        if upper in seen:
            continue
        seen.add(upper)
        ordered.append(upper)
    return ordered


def _markers_from_rankings(
    rankings: Mapping[str, Sequence[str]] | np.ndarray | Mapping[str, Any],
    top_n_markers: int,
) -> dict[str, list[str]]:
    if isinstance(rankings, np.ndarray) and rankings.dtype.names:
        groups = rankings.dtype.names
        derived: dict[str, list[str]] = {}
        for group in groups:
            derived[str(group)] = _unique_ordered(rankings[group][:top_n_markers])
        return derived

    if isinstance(rankings, Mapping):
        derived: dict[str, list[str]] = {}
        for cluster, values in rankings.items():
            derived[str(cluster)] = _unique_ordered(values[:top_n_markers])
        return derived
    raise ValueError("Unsupported structure for rank_genes_groups['names'].")


def _extract_markers(adata: AnnData, cluster_key: str, top_n_markers: int) -> dict[str, list[str]]:
    rankings = adata.uns.get("rank_genes_groups")
    if rankings is None or "names" not in rankings:
        raise ValueError(
            "AnnData object is missing `rank_genes_groups` results for cluster key "
            f"'{cluster_key}'. Run scanpy.tl.rank_genes_groups beforehand or allow "
            "annotate_anndata to compute them."
        )
    return _markers_from_rankings(rankings["names"], top_n_markers)


def _build_cluster_index(cluster_series: pd.Series) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for obs_idx, cluster_id in cluster_series.items():
        mapping.setdefault(str(cluster_id), []).append(obs_idx)
    return mapping


def _report_to_dataframe(report: DatasetReport | dict[str, Any]) -> pd.DataFrame:
    dataset = report if isinstance(report, DatasetReport) else DatasetReport.model_validate(report)
    records: list[dict[str, Any]] = []
    for cluster in dataset.clusters:
        annotation = cluster.annotation
        records.append(
            {
                "cluster_id": cluster.cluster_id,
                "primary_label": annotation.get("primary_label"),
                "confidence": cluster.confidence,
                "status": cluster.status,
                "ontology_id": annotation.get("ontology_id"),
                "rationale": annotation.get("rationale"),
                "warnings": "; ".join(cluster.warnings or []),
            }
        )
    return pd.DataFrame.from_records(records)


def report_to_dataframe(report: DatasetReport | dict[str, Any]) -> pd.DataFrame:
    """Return a dataframe summarising cluster annotations."""

    return _report_to_dataframe(report)


def _resolve_marker_db_path(marker_db_path: Path | None) -> Path:
    if marker_db_path is not None:
        candidate = Path(marker_db_path).expanduser()
        if candidate.exists():
            return candidate

    env_path = os.environ.get(MARKER_DB_ENV_VAR)
    if env_path:
        candidate = Path(env_path).expanduser()
        if candidate.exists():
            return candidate

    settings = get_settings()
    data_dir = Path(settings.data_dir)
    default_path = data_dir / "marker_db.parquet"
    if default_path.exists():
        return default_path

    ensured_dir = assets.ensure_marker_database(target_dir=data_dir)
    candidate = ensured_dir / "marker_db.parquet"
    if candidate.exists():
        return candidate
    fallback = assets.ensure_marker_database() / "marker_db.parquet"
    if fallback.exists():
        return fallback
    raise FileNotFoundError("marker_db.parquet could not be located or materialised.")


_MARKER_CACHE: dict[Path, tuple[float, pd.DataFrame]] = {}


def _load_marker_db(
    marker_db: pd.DataFrame | None,
    marker_db_path: str | Path | None,
    *,
    cache: bool = True,
) -> pd.DataFrame:
    if marker_db is not None:
        missing = [column for column in MARKER_DB_COLUMNS if column not in marker_db.columns]
        if missing:
            raise ValueError(f"Marker DB missing expected columns: {missing}")
        return marker_db[MARKER_DB_COLUMNS].copy()

    resolved = _resolve_marker_db_path(
        Path(marker_db_path) if marker_db_path is not None else None
    )
    resolved = resolved.resolve()

    if cache:
        cached = _MARKER_CACHE.get(resolved)
        mtime = resolved.stat().st_mtime
        if cached and cached[0] == mtime:
            return cached[1].copy()

    df = pd.read_parquet(resolved)
    trimmed = df[MARKER_DB_COLUMNS].copy()
    if cache:
        _MARKER_CACHE[resolved] = (resolved.stat().st_mtime, trimmed)
    return trimmed.copy()


def _guardrail_settings(guardrails: GuardrailConfig | None) -> tuple[Settings | None, int]:
    if guardrails is None:
        settings = get_settings()
        return None, max(1, settings.validation_min_marker_overlap)

    update = guardrails.as_update()
    if not update:
        settings = get_settings()
        return None, max(1, settings.validation_min_marker_overlap)

    base = get_settings()
    override = base.model_copy(update=update)
    min_support = max(1, override.validation_min_marker_overlap)
    return override, min_support


def _build_cache_payload(
    cluster: dict[str, Any],
    dataset_context: dict[str, Any],
    *,
    mode: str,
) -> dict[str, Any]:
    return {
        "version": CACHE_VERSION,
        "cluster_id": str(cluster.get("cluster_id")),
        "markers": cluster.get("markers"),
        "context": dataset_context,
        "mode": mode,
    }


async def _annotate_batches_async(
    annotator: Annotator,
    clusters_payload: list[dict[str, Any]],
    dataset_context: dict[str, Any],
    *,
    batch_options: BatchOptions,
    cache: AnnotationCacheProtocol | None,
    request_id: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    total = len(clusters_payload)
    if total == 0:
        return [], {"total_clusters": 0, "cache_hits": 0, "llm_batches": 0}

    options = batch_options.normalized(total)
    semaphore = asyncio.Semaphore(options.concurrency)
    results: dict[str, dict[str, Any]] = {}
    cache_hits = 0
    llm_batches = 0

    task_group: list[asyncio.Task[None]] = []

    async def process_chunk(index: int, chunk: list[dict[str, Any]]) -> None:
        nonlocal cache_hits, llm_batches
        async with semaphore:
            pending: list[dict[str, Any]] = []
            payloads: list[dict[str, Any]] = []
            if cache is not None:
                for cluster in chunk:
                    payload = _build_cache_payload(
                        cluster,
                        dataset_context,
                        mode=annotator.llm_mode,
                    )
                    cached = await cache.get(payload)
                    if cached:
                        cache_hits += 1
                        results[str(cluster["cluster_id"])] = dict(cached)
                    else:
                        pending.append(cluster)
                        payloads.append(payload)
            else:
                pending = chunk

            if not pending:
                return

            llm_batches += 1
            _log(
                "scanpy.annotate.batch",
                batch_index=index,
                batch_size=len(pending),
                concurrency=options.concurrency,
                request_id=request_id,
            )

            start = time.perf_counter()

            def compute() -> dict[str, Any]:
                return annotator.annotate_batch(pending, dataset_context)

            raw = await _run_in_thread(compute)
            duration = time.perf_counter() - start
            if _BATCH_DURATION is not None:
                _BATCH_DURATION.labels(mode=annotator.llm_mode).observe(duration)
                _BATCH_COUNTER.labels(mode=annotator.llm_mode).inc()
                _CLUSTER_COUNTER.labels(mode=annotator.llm_mode).inc(len(pending))
            for position, cluster in enumerate(pending):
                cluster_id = str(cluster["cluster_id"])
                result = dict(raw.get(cluster_id) or {})
                results[cluster_id] = result
                if cache is not None:
                    payload = payloads[position]
                    await cache.set(payload, result)

    for idx in range(0, total, options.size):
        chunk = clusters_payload[idx : idx + options.size]
        task_group.append(asyncio.create_task(process_chunk(idx // options.size, chunk)))

    if task_group:
        await asyncio.gather(*task_group)

    annotations: list[dict[str, Any]] = []
    for cluster in clusters_payload:
        cluster_id = str(cluster["cluster_id"])
        cluster_result = dict(results.get(cluster_id) or {})
        cluster_result.setdefault("primary_label", "Unknown or Novel")
        cluster_result.setdefault("confidence", "Unknown")
        cluster_result.setdefault("rationale", "")
        annotations.append(
            {
                **cluster_result,
                "cluster_id": cluster_id,
                "markers": cluster.get("markers", []),
            }
        )

    stats = {
        "total_clusters": total,
        "cache_hits": cache_hits,
        "llm_batches": llm_batches,
        "batch_size": options.size,
        "concurrency": options.concurrency,
        "llm_mode": annotator.llm_mode,
    }
    return annotations, stats


async def _run_annotation_workflow(
    clusters_payload: list[dict[str, Any]],
    *,
    annotator: Annotator,
    dataset_context: dict[str, Any],
    marker_df: pd.DataFrame,
    guardrails: GuardrailConfig | None,
    batch_options: BatchOptions,
    cache: AnnotationCacheProtocol | None,
    request_id: str,
) -> tuple[list[dict[str, Any]], DatasetReport, dict[str, Any]]:
    annotations, stats = await _annotate_batches_async(
        annotator,
        clusters_payload,
        dataset_context,
        batch_options=batch_options,
        cache=cache,
        request_id=request_id,
    )

    settings_override, min_support = _guardrail_settings(guardrails)
    crosschecked = await _run_in_thread(
        crosscheck_batch,
        annotations,
        marker_df,
        species=dataset_context.get("species"),
        tissue=dataset_context.get("tissue"),
        min_support=min_support,
    )
    report_model: DatasetReport = build_structured_report(
        annotations,
        crosschecked,
        settings_override=settings_override,
    )
    stats.update(
        {
            "request_id": request_id,
            "min_marker_overlap": min_support,
            "guardrail_overrides": guardrails.as_update() if guardrails else {},
            "summary": report_model.summary.model_dump(),
        }
    )
    return annotations, report_model, stats


def _apply_annotations_to_obs(
    adata: AnnData,
    cluster_key: str,
    result_prefix: str,
    report: DatasetReport,
) -> None:
    label_col = f"{result_prefix}_label"
    confidence_col = f"{result_prefix}_confidence"
    status_col = f"{result_prefix}_status"
    rationale_col = f"{result_prefix}_rationale"
    ontology_col = f"{result_prefix}_ontology_id"
    proposed_col = f"{result_prefix}_proposed_label"
    canonical_col = f"{result_prefix}_canonical_markers"
    mapping_col = f"{result_prefix}_mapping_notes"

    clusters_index = _build_cluster_index(adata.obs[cluster_key].astype(str))
    for column in [
        label_col,
        proposed_col,
        confidence_col,
        status_col,
        rationale_col,
        ontology_col,
        canonical_col,
        mapping_col,
    ]:
        adata.obs[column] = pd.Series(
            ["" for _ in range(adata.n_obs)],
            index=adata.obs.index,
            dtype="object",
        )

    for cluster_report in report.clusters:
        cluster_id = cluster_report.cluster_id
        annotation = cluster_report.annotation
        metadata = annotation.get("metadata") or {}
        mapping_notes = (
            metadata.get("mapping_notes")
            or annotation.get("mapping_notes")
            or []
        )
        canonical_markers = metadata.get("canonical_markers") or annotation.get("markers") or []
        mapping_display = " | ".join(
            f"{note.get('source', '?')}→{note.get('target') or 'unmapped'}" for note in mapping_notes
        )
        obs_indices = clusters_index.get(cluster_id, [])
        for obs_idx in obs_indices:
            adata.obs.at[obs_idx, label_col] = annotation.get("primary_label") or ""
            adata.obs.at[obs_idx, proposed_col] = annotation.get("proposed_label") or ""
            adata.obs.at[obs_idx, confidence_col] = cluster_report.confidence or ""
            adata.obs.at[obs_idx, status_col] = cluster_report.status or ""
            adata.obs.at[obs_idx, rationale_col] = annotation.get("rationale") or ""
            adata.obs.at[obs_idx, ontology_col] = annotation.get("ontology_id") or ""
            adata.obs.at[obs_idx, canonical_col] = ", ".join(canonical_markers)
            adata.obs.at[obs_idx, mapping_col] = mapping_display


async def annotate_anndata_async(
    adata: AnnData,
    cluster_key: str,
    *,
    species: str,
    tissue: str | None = None,
    top_n_markers: int = 5,
    result_prefix: str = "gptca",
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    compute_rankings: bool = True,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: AnnotationCacheProtocol | None = None,
    request_id: str | None = None,
) -> ScanpyAnnotationResult:
    """Asynchronous variant of annotate_anndata."""

    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not found in adata.obs.")

    if compute_rankings:
        _ensure_rankings(
            adata,
            cluster_key,
            top_n_markers=top_n_markers,
        )

    cluster_markers = _extract_markers(adata, cluster_key, top_n_markers)
    clusters_payload = [
        {"cluster_id": cluster_id, "markers": markers}
        for cluster_id, markers in cluster_markers.items()
    ]

    annotator_instance = annotator or Annotator()
    dataset_context: dict[str, str] = {"species": species}
    if tissue:
        dataset_context["tissue"] = tissue

    request = _default_request_id(request_id)
    marker_df = await _run_in_thread(
        _load_marker_db,
        marker_db,
        marker_db_path,
        cache=True,
    )
    annotations, report_model, stats = await _run_annotation_workflow(
        clusters_payload,
        annotator=annotator_instance,
        dataset_context=dataset_context,
        marker_df=marker_df,
        guardrails=guardrails,
        batch_options=batch_options or BatchOptions(),
        cache=annotation_cache,
        request_id=request,
    )

    _apply_annotations_to_obs(adata, cluster_key, result_prefix, report_model)
    return ScanpyAnnotationResult(adata=adata, report=report_model, annotations=annotations, stats=stats)


def annotate_anndata(
    adata: AnnData,
    cluster_key: str,
    *,
    species: str,
    tissue: str | None = None,
    top_n_markers: int = 5,
    result_prefix: str = "gptca",
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    compute_rankings: bool = True,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: AnnotationCacheProtocol | None = None,
    request_id: str | None = None,
) -> ScanpyAnnotationResult:
    """Annotate clusters in an AnnData object using GPT Cell Annotator."""

    try:
        loop = asyncio.get_running_loop()
        if loop.is_running():
            raise RuntimeError(
                "annotate_anndata cannot be called from an async context. "
                "Use `await annotate_anndata_async(...)` instead."
            )
    except RuntimeError:
        loop = None

    return asyncio.run(
        annotate_anndata_async(
            adata,
            cluster_key,
            species=species,
            tissue=tissue,
            top_n_markers=top_n_markers,
            result_prefix=result_prefix,
            marker_db=marker_db,
            marker_db_path=marker_db_path,
            annotator=annotator,
            compute_rankings=compute_rankings,
            batch_options=batch_options,
            guardrails=guardrails,
            annotation_cache=annotation_cache,
            request_id=request_id,
        )
    )


async def annotate_from_markers_async(
    cluster_markers: Mapping[str, Sequence[str]],
    *,
    species: str,
    tissue: str | None = None,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: AnnotationCacheProtocol | None = None,
    request_id: str | None = None,
) -> MarkerAnnotationResult:
    """Annotate clusters from a pre-computed marker mapping."""

    clusters_payload = [
        {"cluster_id": str(cluster_id), "markers": _unique_ordered(markers)}
        for cluster_id, markers in cluster_markers.items()
    ]
    annotator_instance = annotator or Annotator()
    dataset_context: dict[str, str] = {"species": species}
    if tissue:
        dataset_context["tissue"] = tissue
    marker_df = await _run_in_thread(
        _load_marker_db,
        marker_db,
        marker_db_path,
        cache=True,
    )
    annotations, report_model, stats = await _run_annotation_workflow(
        clusters_payload,
        annotator=annotator_instance,
        dataset_context=dataset_context,
        marker_df=marker_df,
        guardrails=guardrails,
        batch_options=batch_options or BatchOptions(),
        cache=annotation_cache,
        request_id=_default_request_id(request_id),
    )
    return MarkerAnnotationResult(annotations=annotations, report=report_model, stats=stats)


def annotate_from_markers(
    cluster_markers: Mapping[str, Sequence[str]],
    *,
    species: str,
    tissue: str | None = None,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: AnnotationCacheProtocol | None = None,
    request_id: str | None = None,
) -> MarkerAnnotationResult:
    return asyncio.run(
        annotate_from_markers_async(
            cluster_markers,
            species=species,
            tissue=tissue,
            marker_db=marker_db,
            marker_db_path=marker_db_path,
            annotator=annotator,
            batch_options=batch_options,
            guardrails=guardrails,
            annotation_cache=annotation_cache,
            request_id=request_id,
        )
    )


def annotate_rank_genes(
    rank_genes_result: Mapping[str, Sequence[str]] | Mapping[str, Any],
    *,
    species: str,
    tissue: str | None = None,
    top_n_markers: int = 5,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: AnnotationCacheProtocol | None = None,
    request_id: str | None = None,
) -> MarkerAnnotationResult:
    cluster_markers = _markers_from_rankings(rank_genes_result, top_n_markers)
    return annotate_from_markers(
        cluster_markers,
        species=species,
        tissue=tissue,
        marker_db=marker_db,
        marker_db_path=marker_db_path,
        annotator=annotator,
        batch_options=batch_options,
        guardrails=guardrails,
        annotation_cache=annotation_cache,
        request_id=request_id,
    )


def _cluster_labels_from_obs(
    adata: AnnData,
    cluster_key: str,
    label_column: str,
) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for obs_idx, cluster_id in adata.obs[cluster_key].astype(str).items():
        label = adata.obs.at[obs_idx, label_column]
        if isinstance(label, str) and label.strip():
            mapping.setdefault(cluster_id, label.strip())
    return mapping


def validate_anndata(
    adata: AnnData,
    cluster_key: str,
    *,
    species: str,
    label_column: str,
    ontology_column: str | None = None,
    tissue: str | None = None,
    top_n_markers: int = 5,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    compute_rankings: bool = True,
    guardrails: GuardrailConfig | None = None,
) -> DatasetReport:
    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not found in adata.obs.")
    if label_column not in adata.obs:
        raise KeyError(f"Label column '{label_column}' not found in adata.obs.")
    if ontology_column and ontology_column not in adata.obs:
        raise KeyError(f"Ontology column '{ontology_column}' not found in adata.obs.")

    if compute_rankings:
        _ensure_rankings(adata, cluster_key, top_n_markers=top_n_markers)

    cluster_markers = _extract_markers(adata, cluster_key, top_n_markers)
    label_mapping = _cluster_labels_from_obs(adata, cluster_key, label_column)
    ontology_mapping = (
        _cluster_labels_from_obs(adata, cluster_key, ontology_column) if ontology_column else {}
    )

    annotations: list[dict[str, Any]] = []
    for cluster_id, markers in cluster_markers.items():
        annotation = {
            "cluster_id": cluster_id,
            "markers": markers,
            "primary_label": label_mapping.get(cluster_id, "Unknown or Novel"),
        }
        if ontology_column:
            annotation["ontology_id"] = ontology_mapping.get(cluster_id)
        annotations.append(annotation)

    marker_df = _load_marker_db(marker_db, marker_db_path, cache=True)
    settings_override, min_support = _guardrail_settings(guardrails)
    crosschecked = crosscheck_batch(
        annotations,
        marker_df,
        species=species,
        tissue=tissue,
        min_support=min_support,
    )
    return build_structured_report(
        annotations,
        crosschecked,
        settings_override=settings_override,
    )


def _read_adata(path: Path) -> AnnData:
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        return ad.read_h5ad(path)
    if suffix == ".loom":
        if sc is None:
            raise ImportError("scanpy is required to read Loom files.")
        return sc.read_loom(path, obsm_keys=None)
    raise ValueError(f"Unsupported AnnData format: {path.suffix}")


def _write_report(report: DatasetReport, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report.model_dump(), indent=2), encoding="utf-8")


def _apply_preset(args: argparse.Namespace) -> None:
    if not args.preset:
        return
    preset = SPECIES_PRESETS.get(args.preset)
    if not preset:
        raise SystemExit(f"Unknown preset '{args.preset}'. Available presets: {', '.join(SPECIES_PRESETS)}")
    if not args.species:
        args.species = preset.get("species")
    if not args.tissue and preset.get("tissue"):
        args.tissue = preset["tissue"]


def _annotation_cache_from_args(args: argparse.Namespace) -> AnnotationCacheProtocol | None:
    if args.cache_dir:
        return DiskAnnotationCache(Path(args.cache_dir).expanduser())

    env_dir = os.environ.get(CACHE_DIR_ENV_VAR)
    if env_dir:
        return DiskAnnotationCache(Path(env_dir).expanduser())
    return None


def _build_batch_options(args: argparse.Namespace) -> BatchOptions:
    return BatchOptions(size=args.batch_size, concurrency=args.concurrency)


def cmd_annotate(args: argparse.Namespace) -> int:
    _apply_preset(args)
    cache = _annotation_cache_from_args(args)
    annotator = Annotator(force_mock=args.offline) if args.offline else Annotator()
    adata = _read_adata(args.input)
    result = annotate_anndata(
        adata,
        args.cluster_key,
        species=args.species,
        tissue=args.tissue,
        top_n_markers=args.top_n_markers,
        result_prefix=args.prefix,
        marker_db_path=args.marker_db,
        annotator=annotator,
        compute_rankings=not args.skip_recompute_markers,
        batch_options=_build_batch_options(args),
        guardrails=GuardrailConfig(
            min_marker_overlap=args.guardrail_min_overlap,
            force_unknown_on_fail=args.guardrail_force_unknown,
        )
        if args.guardrail_min_overlap is not None or args.guardrail_force_unknown is not None
        else None,
        annotation_cache=cache,
        request_id=args.request_id,
    )

    output_path = args.output or args.input
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.adata.write(output_path)

    if args.summary_csv:
        df = report_to_dataframe(result.report)
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.summary_csv, index=False)
    if args.summary_json:
        _write_report(result.report, args.summary_json)
    if args.stats_json:
        args.stats_json.parent.mkdir(parents=True, exist_ok=True)
        args.stats_json.write_text(json.dumps(result.stats, indent=2), encoding="utf-8")
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    _apply_preset(args)
    adata = _read_adata(args.input)
    report = validate_anndata(
        adata,
        args.cluster_key,
        species=args.species,
        tissue=args.tissue,
        label_column=args.label_column,
        ontology_column=args.ontology_column,
        top_n_markers=args.top_n_markers,
        marker_db_path=args.marker_db,
        compute_rankings=not args.skip_recompute_markers,
        guardrails=GuardrailConfig(
            min_marker_overlap=args.guardrail_min_overlap,
            force_unknown_on_fail=args.guardrail_force_unknown,
        )
        if args.guardrail_min_overlap is not None or args.guardrail_force_unknown is not None
        else None,
    )
    if args.summary_csv:
        df = report_to_dataframe(report)
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.summary_csv, index=False)
    if args.summary_json:
        _write_report(report, args.summary_json)
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m gpt_cell_annotator.scanpy",
        description="Scanpy integration tools for GPT Cell Annotator.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate clusters in an AnnData file (h5ad or loom).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    annotate_parser.add_argument("input", type=Path, help="Path to the input AnnData file.")
    annotate_parser.add_argument(
        "--cluster-key",
        required=True,
        help="Column in adata.obs that identifies cluster assignments.",
    )
    annotate_parser.add_argument("--species", help="Species context (overrides preset).")
    annotate_parser.add_argument("--preset", choices=sorted(SPECIES_PRESETS), help="Use preset species/tissue.")
    annotate_parser.add_argument("--tissue", help="Optional tissue context.")
    annotate_parser.add_argument(
        "--top-n-markers",
        type=int,
        default=5,
        help="Number of top markers per cluster to send for annotation.",
    )
    annotate_parser.add_argument(
        "--prefix",
        default="gptca",
        help="Prefix for columns written to adata.obs.",
    )
    annotate_parser.add_argument(
        "--output",
        type=Path,
        help="Destination AnnData path. Defaults to overwriting input.",
    )
    annotate_parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Optional path to write a per-cluster summary CSV.",
    )
    annotate_parser.add_argument(
        "--summary-json",
        type=Path,
        help="Optional path to write the full annotation report JSON.",
    )
    annotate_parser.add_argument(
        "--stats-json",
        type=Path,
        help="Optional path to write batching and telemetry statistics.",
    )
    annotate_parser.add_argument(
        "--marker-db",
        type=Path,
        help="Path to marker_db.parquet. Defaults to env var or cached assets.",
    )
    annotate_parser.add_argument(
        "--skip-recompute-markers",
        action="store_true",
        help="Assume rank_genes_groups already computed; do not recompute.",
    )
    annotate_parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Number of clusters per annotator batch.",
    )
    annotate_parser.add_argument(
        "--concurrency",
        type=int,
        default=1,
        help="Maximum concurrent batches (uses asyncio threads).",
    )
    annotate_parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Directory for disk cache of annotation responses.",
    )
    annotate_parser.add_argument(
        "--offline",
        action="store_true",
        help="Force heuristic mock mode and avoid network access.",
    )
    annotate_parser.add_argument(
        "--guardrail-min-overlap",
        type=int,
        help="Override validation_min_marker_overlap for this run.",
    )
    annotate_parser.add_argument(
        "--guardrail-force-unknown",
        action="store_true",
        help="Force downgrade to Unknown when guardrails fail.",
    )
    annotate_parser.add_argument(
        "--request-id",
        help="Attach a custom request identifier for telemetry.",
    )
    annotate_parser.set_defaults(func=cmd_annotate)

    validate_parser = subparsers.add_parser(
        "validate",
        help="Run guardrail validation on existing annotations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    validate_parser.add_argument("input", type=Path, help="Path to the input AnnData file.")
    validate_parser.add_argument(
        "--cluster-key",
        required=True,
        help="Column in adata.obs that identifies cluster assignments.",
    )
    validate_parser.add_argument("--species", required=True, help="Species context.")
    validate_parser.add_argument("--preset", choices=sorted(SPECIES_PRESETS), help="Use preset species/tissue.")
    validate_parser.add_argument("--tissue", help="Optional tissue context.")
    validate_parser.add_argument(
        "--label-column",
        required=True,
        help="Existing annotation column to validate.",
    )
    validate_parser.add_argument(
        "--ontology-column",
        help="Optional ontology column to validate.",
    )
    validate_parser.add_argument(
        "--top-n-markers",
        type=int,
        default=5,
        help="Number of top markers per cluster for guardrail checks.",
    )
    validate_parser.add_argument(
        "--marker-db",
        type=Path,
        help="Path to marker_db.parquet. Defaults to env var or cached assets.",
    )
    validate_parser.add_argument(
        "--skip-recompute-markers",
        action="store_true",
        help="Assume rank_genes_groups already computed; do not recompute.",
    )
    validate_parser.add_argument(
        "--guardrail-min-overlap",
        type=int,
        help="Override validation_min_marker_overlap for this run.",
    )
    validate_parser.add_argument(
        "--guardrail-force-unknown",
        action="store_true",
        help="Force downgrade to Unknown when guardrails fail.",
    )
    validate_parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Optional path to write per-cluster validation summary CSV.",
    )
    validate_parser.add_argument(
        "--summary-json",
        type=Path,
        help="Optional path to write validation report JSON.",
    )
    validate_parser.set_defaults(func=cmd_validate)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    return args.func(args)


__all__ = [
    "BatchOptions",
    "DiskAnnotationCache",
    "GuardrailConfig",
    "MarkerAnnotationResult",
    "ScanpyAnnotationResult",
    "annotate_anndata",
    "annotate_anndata_async",
    "annotate_from_markers",
    "annotate_from_markers_async",
    "annotate_rank_genes",
    "report_to_dataframe",
    "validate_anndata",
    "main",
    "MARKER_DB_COLUMNS",
]
