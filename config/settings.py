"""Application configuration via environment variables."""

from __future__ import annotations

from functools import lru_cache

from pydantic import Field
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Central configuration for CellAnnot-GPT services."""

    openai_api_key: str = Field(default="", alias="OPENAI_API_KEY")
    openai_model: str = Field(default="gpt-4o-mini")
    openai_temperature: float = Field(default=0.2)
    openai_max_tokens: int = Field(default=800)
    openai_requests_per_minute: int = Field(default=20)
    openai_retry_attempts: int = Field(default=3)
    openai_retry_backoff_seconds: float = Field(default=1.5)

    environment: str = Field(default="development")
    data_dir: str = Field(default="data/processed")
    redis_url: str = Field(default="", alias="REDIS_URL")

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"
        populate_by_name = True
        extra = "ignore"


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return cached settings instance."""

    return Settings()


__all__ = ["Settings", "get_settings"]
