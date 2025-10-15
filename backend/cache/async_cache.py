"""Async cache helper for annotator integration."""

from __future__ import annotations

import asyncio
from collections.abc import Awaitable, Callable
from typing import Any

from backend.cache.redis_cache import RedisAnnotationCache


async def get_cached_or_run(
    cache: RedisAnnotationCache | None,
    payload: dict[str, Any],
    compute: Callable[[], Awaitable[Any]] | Callable[[], Any],
) -> tuple[Any, bool]:
    if cache is None:
        return await _evaluate(compute), False

    cached = await cache.get(payload)
    if cached is not None:
        return cached, True

    result = await _evaluate(compute)
    await cache.set(payload, result)
    return result, False


async def _evaluate(func: Callable[[], Awaitable[Any]] | Callable[[], Any]) -> Any:
    if asyncio.iscoroutinefunction(func):  # type: ignore[arg-type]
        return await func()  # type: ignore[misc]
    result = func()
    if asyncio.iscoroutine(result):
        return await result  # type: ignore[misc]
    return await asyncio.to_thread(lambda: result)
