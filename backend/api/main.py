"""FastAPI application exposing CellAnnot-GPT endpoints."""

from __future__ import annotations

import logging
from typing import Dict, Iterable, Optional

from fastapi import Depends, FastAPI, HTTPException, status

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
import pandas as pd

logger = logging.getLogger("cellannot.api")


def get_annotator() -> Annotator:
    """Dependency provider for the annotator."""

    return Annotator()


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
async def health() -> HealthResponse:
    return HealthResponse(status="ok")


@app.post("/annotate_cluster", response_model=AnnotationResponse)
async def annotate_cluster(
    payload: AnnotateClusterRequest,
    annotator: Annotator = Depends(get_annotator),
    marker_db: pd.DataFrame = Depends(get_marker_db),
) -> AnnotationResponse:
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

    crosschecked = crosscheck_batch(
        [annotation_record],
        marker_db,
        species=(payload.dataset_context.species if payload.dataset_context else None),
        tissue=(payload.dataset_context.tissue if payload.dataset_context else None),
    )
    report = build_structured_report([annotation_record], crosschecked)
    return AnnotationResponse(result=report)


@app.post("/annotate_batch", response_model=BatchAnnotationResponse)
async def annotate_batch(
    payload: AnnotateBatchRequest,
    annotator: Annotator = Depends(get_annotator),
    marker_db: pd.DataFrame = Depends(get_marker_db),
) -> BatchAnnotationResponse:
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
    report = build_structured_report(annotations, crosschecked)
    return BatchAnnotationResponse(result=report)


__all__ = ["app"]
