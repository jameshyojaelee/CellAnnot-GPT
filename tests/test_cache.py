from __future__ import annotations

import json

import pytest

from backend.cache.async_cache import get_cached_or_run
from backend.cache.redis_cache import CacheConfig, RedisAnnotationCache

try:
    import fakeredis.aioredis as fakeredis
except ImportError:  # pragma: no cover - optional dependency
    fakeredis = None  # type: ignore


@pytest.mark.skipif(fakeredis is None, reason="fakeredis not installed")
@pytest.mark.asyncio
async def test_cache_roundtrip(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_client = fakeredis.FakeRedis(decode_responses=True)

    monkeypatch.setattr(
        "backend.cache.redis_cache.redis.from_url",
        lambda *args, **kwargs: fake_client,
    )
    cache = RedisAnnotationCache(CacheConfig(url="redis://localhost"))

    payload = {"type": "cluster", "markers": ["A"], "context": {}, "llm_mode": "mock"}
    value = {"summary": {"total_clusters": 1}}

    await cache.set(payload, value)
    cached = await cache.get(payload)
    assert cached == value


@pytest.mark.asyncio
async def test_get_cached_or_run_without_cache_async():
    async def compute():
        return {"value": 42}

    result, from_cache = await get_cached_or_run(None, {"key": "a"}, compute)
    assert not from_cache
    assert result == {"value": 42}


@pytest.mark.asyncio
async def test_get_cached_or_run_with_cache_and_sync_compute():
    class MemoryCache:
        def __init__(self) -> None:
            self.store = {}

        async def get(self, payload):
            return self.store.get(json.dumps(payload, sort_keys=True))

        async def set(self, payload, value):
            self.store[json.dumps(payload, sort_keys=True)] = value

    cache = MemoryCache()
    calls = {"count": 0}

    def compute():
        calls["count"] += 1
        return {"result": "fresh"}

    payload = {"type": "cluster", "id": "1"}

    first, from_cache_first = await get_cached_or_run(cache, payload, compute)
    assert not from_cache_first
    assert first == {"result": "fresh"}
    assert calls["count"] == 1

    second, from_cache_second = await get_cached_or_run(cache, payload, compute)
    assert from_cache_second
    assert second == {"result": "fresh"}
    assert calls["count"] == 1
