"""Async cache helper for annotator integration."""

from __future__ import annotations

import asyncio
from typing import Any, Dict, Optional

from backend.cache.redis_cache import RedisAnnotationCache


async def get_cached_or_run(
    cache: Optional[RedisAnnotationCache],
    payload: Dict[str, Any],
    compute: callable,
):  # type: ignore[not-callable]
    if cache is None:
        return await asyncio.to_thread(compute), False

    cached = await cache.get(payload)
    if cached is not None:
        return cached, True

    result = await asyncio.to_thread(compute)
    await cache.set(payload, result)
    return result, False
