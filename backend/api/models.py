"""Request/response models for the GPT Cell Annotator API."""

from __future__ import annotations

from pydantic import BaseModel, Field


class ClusterStats(BaseModel):
    avg_log_fc: float | None = Field(None, alias="avg_logFC")
    pct_expressed: float | None


class ClusterPayload(BaseModel):
    cluster_id: str
    markers: list[str] = Field(default_factory=list)
    stats: ClusterStats | None = None
    notes: str | None = None


class DatasetContext(BaseModel):
    species: str | None = None
    tissue: str | None = None
    condition: str | None = None


class AnnotateClusterRequest(BaseModel):
    cluster: ClusterPayload
    dataset_context: DatasetContext | None = None
    return_validated: bool = False


class AnnotateBatchRequest(BaseModel):
    clusters: list[ClusterPayload]
    dataset_context: DatasetContext | None = None


class AnnotationResponse(BaseModel):
    result: dict


class BatchAnnotationResponse(BaseModel):
    result: dict


class HealthResponse(BaseModel):
    status: str = "ok"
    llm_mode: str = "unknown"
    cache_enabled: bool = False


__all__ = [
    "AnnotateBatchRequest",
    "AnnotateClusterRequest",
    "AnnotationResponse",
    "BatchAnnotationResponse",
    "ClusterPayload",
    "ClusterStats",
    "DatasetContext",
    "HealthResponse",
]
