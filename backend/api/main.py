"""FastAPI application exposing CellAnnot-GPT endpoints."""

from __future__ import annotations

import logging
from typing import Dict, Iterable, Optional

from fastapi import Depends, FastAPI, HTTPException, status

from backend.api.logging_config import configure_logging
from backend.api.middleware import logging_middleware
from backend.api.models import (
    AnnotateBatchRequest,
    AnnotateClusterRequest,
    AnnotationResponse,
    BatchAnnotationResponse,
    HealthResponse,
)
from backend.llm.annotator import Annotator, AnnotationError
from backend.validation.crosscheck import crosscheck_batch
from backend.validation.report import build_structured_report
from backend.data_ingest.marker_loader import NORMALIZED_COLUMNS
from backend.cache.redis_cache import RedisAnnotationCache, create_cache_from_env
from config.settings import get_settings
import pandas as pd

logger = logging.getLogger("cellannot.api")
settings = get_settings()
_CACHE = create_cache_from_env(settings.redis_url)

configure_logging(settings.log_level if hasattr(settings, "log_level") else "INFO")


def get_annotator() -> Annotator:
    """Dependency provider for the annotator."""

    return Annotator()


def get_cache() -> Optional[RedisAnnotationCache]:
    return _CACHE


def get_marker_db() -> pd.DataFrame:
    """Load marker database from default artifact location."""

    try:
        df = pd.read_parquet("data/processed/marker_db.parquet")
        return df[NORMALIZED_COLUMNS]
    except FileNotFoundError as exc:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Marker database not built. Run scripts/build_marker_db.py first.",
        ) from exc


app = FastAPI(title="CellAnnot-GPT API")
app.middleware("http")(logging_middleware)


@app.get("/health", response_model=HealthResponse)
async def health(
    annotator: Annotator = Depends(get_annotator),
    cache: Optional[RedisAnnotationCache] = Depends(get_cache),
) -> HealthResponse:
    return HealthResponse(status="ok", llm_mode=annotator.llm_mode, cache_enabled=cache is not None)


@app.post("/annotate_cluster", response_model=AnnotationResponse)
async def annotate_cluster(
    payload: AnnotateClusterRequest,
    annotator: Annotator = Depends(get_annotator),
    marker_db: pd.DataFrame = Depends(get_marker_db),
    cache: Optional[RedisAnnotationCache] = Depends(get_cache),
) -> AnnotationResponse:
    logger.info(
        "Annotating cluster %s using llm_mode=%s",
        payload.cluster.cluster_id,
        annotator.llm_mode,
    )
    context_dict = payload.dataset_context.model_dump(exclude_none=True) if payload.dataset_context else {}
    cache_payload = {
        "type": "cluster",
        "markers": sorted(payload.cluster.markers),
        "context": context_dict,
        "llm_mode": annotator.llm_mode,
        "validated": payload.return_validated,
    }
    if cache is not None and not payload.return_validated:
        try:
            cached = await cache.get(cache_payload)
        except Exception as exc:  # pragma: no cover - redis availability issue
            logger.warning("Redis get failed: %s", exc)
        else:
            if cached is not None:
                logger.debug("Cache hit for cluster %s", payload.cluster.cluster_id)
                return AnnotationResponse(result=cached)
    try:
        result = annotator.annotate_cluster(
            payload.cluster.model_dump(),
            payload.dataset_context.model_dump(exclude_none=True) if payload.dataset_context else None,
        )
    except AnnotationError as exc:
        logger.exception("Annotation failed")
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    annotation_record = {
        **result,
        "cluster_id": payload.cluster.cluster_id,
        "markers": payload.cluster.markers,
    }

    if payload.return_validated:
        crosschecked = crosscheck_batch(
            [annotation_record],
            marker_db,
            species=(payload.dataset_context.species if payload.dataset_context else None),
            tissue=(payload.dataset_context.tissue if payload.dataset_context else None),
        )
        report_model = build_structured_report([annotation_record], crosschecked)
        report = report_model.model_dump()
        return AnnotationResponse(result=report)

    if cache is not None:
        try:
            await cache.set(cache_payload, annotation_record)
        except Exception as exc:  # pragma: no cover - redis availability issue
            logger.warning("Redis set failed: %s", exc)
    return AnnotationResponse(result=annotation_record)


@app.post("/annotate_batch", response_model=BatchAnnotationResponse)
async def annotate_batch(
    payload: AnnotateBatchRequest,
    annotator: Annotator = Depends(get_annotator),
    marker_db: pd.DataFrame = Depends(get_marker_db),
    cache: Optional[RedisAnnotationCache] = Depends(get_cache),
) -> BatchAnnotationResponse:
    logger.info(
        "Annotating batch of %s clusters using llm_mode=%s",
        len(payload.clusters),
        annotator.llm_mode,
    )
    context_dict = payload.dataset_context.model_dump(exclude_none=True) if payload.dataset_context else {}
    clusters_payload = [
        {"cluster_id": str(cluster.cluster_id), "markers": sorted(cluster.markers)}
        for cluster in payload.clusters
    ]
    cache_payload = {
        "type": "batch",
        "clusters": clusters_payload,
        "context": context_dict,
        "llm_mode": annotator.llm_mode,
        "validated": True,
    }
    if cache is not None:
        try:
            cached = await cache.get(cache_payload)
        except Exception as exc:  # pragma: no cover
            logger.warning("Redis get failed: %s", exc)
        else:
            if cached is not None:
                logger.debug("Cache hit for batch of %s clusters", len(payload.clusters))
                return BatchAnnotationResponse(result=cached)
    try:
        result = annotator.annotate_batch(
            [cluster.model_dump() for cluster in payload.clusters],
            payload.dataset_context.model_dump(exclude_none=True) if payload.dataset_context else None,
        )
    except AnnotationError as exc:
        logger.exception("Batch annotation failed")
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    annotations = []
    for cluster in payload.clusters:
        cluster_result = result.get(str(cluster.cluster_id)) or {}
        annotations.append(
            {
                **cluster_result,
                "cluster_id": cluster.cluster_id,
                "markers": cluster.markers,
            }
        )

    crosschecked = crosscheck_batch(
        annotations,
        marker_db,
        species=(payload.dataset_context.species if payload.dataset_context else None),
        tissue=(payload.dataset_context.tissue if payload.dataset_context else None),
    )
    report_model = build_structured_report(annotations, crosschecked)
    report = report_model.model_dump()
    if cache is not None:
        try:
            await cache.set(cache_payload, report)
        except Exception as exc:  # pragma: no cover
            logger.warning("Redis set failed: %s", exc)
    return BatchAnnotationResponse(result=report)


__all__ = ["app"]
