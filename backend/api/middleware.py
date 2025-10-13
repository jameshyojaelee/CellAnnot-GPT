"""Custom middleware for logging and timing API requests."""

from __future__ import annotations

import logging
import time
from typing import Callable

from fastapi import Request, Response

logger = logging.getLogger("cellannot.api")


async def logging_middleware(request: Request, call_next: Callable[[Request], Response]) -> Response:
    """Log incoming requests with latency in milliseconds."""

    start = time.perf_counter()
    response = await call_next(request)
    duration_ms = (time.perf_counter() - start) * 1000
    logger.info(
        "%s %s %s %.2fms",
        request.method,
        request.url.path,
        response.status_code,
        duration_ms,
    )
    response.headers["X-Process-Time-ms"] = f"{duration_ms:.2f}"
    return response


__all__ = ["logging_middleware"]
