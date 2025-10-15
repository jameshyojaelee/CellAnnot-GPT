"""Redis cache utilities for annotation responses."""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from hashlib import sha256
from typing import Any

try:
    import redis.asyncio as redis
except ImportError:  # pragma: no cover - optional dependency
    redis = None  # type: ignore

logger = logging.getLogger("cellannot.cache")


@dataclass
class CacheConfig:
    url: str
    namespace: str = "cellannot"
    ttl_seconds: int = 60 * 60


class RedisAnnotationCache:
    def __init__(self, config: CacheConfig) -> None:
        if redis is None:
            raise RuntimeError("redis.asyncio is not installed")
        self.config = config
        self._client = redis.from_url(config.url, decode_responses=True)

    def _key(self, payload: dict[str, Any]) -> str:
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        digest = sha256(canonical.encode("utf-8")).hexdigest()
        return f"{self.config.namespace}:{digest}"

    async def get(self, payload: dict[str, Any]) -> dict[str, Any] | None:
        key = self._key(payload)
        raw = await self._client.get(key)
        if raw is None:
            return None
        try:
            decoded = json.loads(raw)
            if isinstance(decoded, dict):
                return decoded
            logger.warning("Unexpected payload type in cache for key %s; evicting", key)
            await self._client.delete(key)
            return None
        except json.JSONDecodeError:
            logger.warning("Invalid JSON in cache for key %s; evicting", key)
            await self._client.delete(key)
            return None

    async def set(self, payload: dict[str, Any], value: dict[str, Any]) -> None:
        key = self._key(payload)
        await self._client.set(key, json.dumps(value), ex=self.config.ttl_seconds)

    async def close(self) -> None:
        await self._client.close()


def create_cache_from_env(url: str | None) -> RedisAnnotationCache | None:
    if not url:
        return None
    if redis is None:
        logger.warning("redis.asyncio not installed; cache disabled")
        return None
    try:
        cache = RedisAnnotationCache(CacheConfig(url=url))
        return cache
    except Exception as exc:  # pragma: no cover - connection failure
        logger.warning("Failed to initialize Redis cache: %s", exc)
        return None


__all__ = ["RedisAnnotationCache", "CacheConfig", "create_cache_from_env"]
