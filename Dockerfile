# syntax=docker/dockerfile:1

FROM python:3.11-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

ENV POETRY_VERSION=1.8.3
RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="/root/.local/bin:$PATH"

WORKDIR /app
COPY pyproject.toml README.md .
COPY poetry.lock poetry.lock
RUN poetry config virtualenvs.in-project true && poetry install --no-root --no-interaction --only main

COPY backend backend
COPY frontend frontend
COPY evaluation evaluation
COPY config config
COPY scripts scripts
COPY schemas schemas
COPY docs docs
COPY prompts.md prompts.md
COPY data data

RUN poetry install --no-interaction

FROM python:3.11-slim AS runtime
WORKDIR /app
COPY --from=builder /app /app
ENV PATH="/app/.venv/bin:$PATH"
ENV PYTHONUNBUFFERED=1
EXPOSE 8000
CMD ["uvicorn", "backend.api.main:app", "--host", "0.0.0.0", "--port", "8000"]
