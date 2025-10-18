# syntax=docker/dockerfile:1

FROM python:3.11-slim AS builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY pyproject.toml poetry.lock README.md ./
COPY backend backend
COPY config config
COPY data data
COPY docs docs
COPY evaluation evaluation
COPY frontend frontend
COPY gpt_cell_annotator gpt_cell_annotator
COPY schemas schemas
COPY scripts scripts
COPY tests tests
RUN pip install --upgrade pip build \
    && python -m build --wheel --sdist

FROM python:3.11-slim AS runtime

ENV PIP_NO_CACHE_DIR=1 \
    PYTHONUNBUFFERED=1 \
    GPT_CELL_ANNOTATOR_HOME=/data \
    GPT_CELL_ANNOTATOR_DATA_DIR=/data/processed

WORKDIR /app
COPY --from=builder /build/dist /tmp/dist
RUN pip install --upgrade pip \
    && python - <<'PY'
import glob
import subprocess
import sys

wheels = sorted(glob.glob("/tmp/dist/gpt_cell_annotator-*.whl"))
if not wheels:
    raise SystemExit("No GPT Cell Annotator wheel found in /tmp/dist")
wheel = wheels[-1]
subprocess.check_call([sys.executable, "-m", "pip", "install", f"{wheel}[api,scanpy]"])
PY
    && rm -rf /tmp/dist

VOLUME ["/data"]
EXPOSE 8000
ENTRYPOINT ["gca"]
CMD ["--help"]
