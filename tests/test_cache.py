from __future__ import annotations

import pytest

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
