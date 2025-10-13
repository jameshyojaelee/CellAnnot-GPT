"""Custom middleware for logging, timing, and tracing HTTP requests."""

from __future__ import annotations

import time
import uuid
from typing import Callable

import structlog
from fastapi import Request, Response

logger = structlog.get_logger("cellannot.api")


async def logging_middleware(request: Request, call_next: Callable[[Request], Response]) -> Response:
    """Attach trace ID and emit structured timing logs."""

    trace_id = request.headers.get("X-Trace-Id", str(uuid.uuid4()))
    start = time.perf_counter()
    structlog.contextvars.clear_contextvars()
    structlog.contextvars.bind_contextvars(trace_id=trace_id, path=request.url.path, method=request.method)

    response = await call_next(request)

    duration_ms = (time.perf_counter() - start) * 1000
    structlog.contextvars.bind_contextvars(status_code=response.status_code, duration_ms=duration_ms)
    logger.info("request.completed")

    response.headers["X-Trace-Id"] = trace_id
    response.headers["X-Process-Time-ms"] = f"{duration_ms:.2f}"
    return response


__all__ = ["logging_middleware"]
