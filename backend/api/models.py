"""Request/response models for the CellAnnot-GPT API."""

from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, Field


class ClusterStats(BaseModel):
    avg_log_fc: Optional[float] = Field(None, alias="avg_logFC")
    pct_expressed: Optional[float]


class ClusterPayload(BaseModel):
    cluster_id: str
    markers: List[str] = Field(default_factory=list)
    stats: Optional[ClusterStats] = None
    notes: Optional[str] = None


class DatasetContext(BaseModel):
    species: Optional[str] = None
    tissue: Optional[str] = None
    condition: Optional[str] = None


class AnnotateClusterRequest(BaseModel):
    cluster: ClusterPayload
    dataset_context: Optional[DatasetContext] = None


class AnnotateBatchRequest(BaseModel):
    clusters: List[ClusterPayload]
    dataset_context: Optional[DatasetContext] = None


class AnnotationResponse(BaseModel):
    result: dict


class BatchAnnotationResponse(BaseModel):
    result: dict


class HealthResponse(BaseModel):
    status: str = "ok"


__all__ = [
    "ClusterPayload",
    "ClusterStats",
    "DatasetContext",
    "AnnotateClusterRequest",
    "AnnotateBatchRequest",
    "AnnotationResponse",
    "BatchAnnotationResponse",
    "HealthResponse",
]
