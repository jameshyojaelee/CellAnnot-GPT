# Deployment Guide

This document describes how to deploy CellAnnot-GPT across development, staging, and production environments, and how to manage runtime secrets securely.

## 1. Runtime Components
- **API service** (`backend.api.main:app` via Uvicorn) – serves annotation endpoints.
- **Redis (optional)** – backing store for caching prompt responses. Included in `docker-compose.yml` but not required for local development.
- **Marker database artifacts** – mounted from `data/processed` so the container can read prebuilt marker knowledge stores.

## 2. Environment Configuration
Define the following environment variables:

| Variable | Description | Required |
| --- | --- | --- |
| `OPENAI_API_KEY` | API key for OpenAI GPT-4o (or compatible model). | Yes (unless running fully offline) |
| `CELLANNOT_DATA_DIR` | Directory containing `marker_db.parquet` / `marker_db.sqlite`. Defaults to `/data`. | Recommended |
| `REDIS_URL` | Redis connection string if caching is enabled. | Optional |
| `ENVIRONMENT` | One of `development`, `staging`, `production` – consumed by `config/settings.py`. | Optional |

### Local `.env`
Create a `.env` file at repo root (not committed) when using Docker Compose:

```env
OPENAI_API_KEY=sk-...
ENVIRONMENT=development
```

## 3. Development Deployment
1. Build + run services with Docker Compose:
   ```bash
   docker compose up --build
   ```
2. The API is available at `http://localhost:8000`. Streamlit can be run separately on the host or added as another service if needed.
3. Marker DB volume `./data/processed:/data` lets you rebuild the knowledge base on the host (`poetry run python scripts/build_marker_db.py`) and reuse it inside the container.

## 4. Staging Deployment
- Build and tag the image:
  ```bash
  docker build -t ghcr.io/<org>/cellannot-gpt:<sha> .
  ```
- Push to your registry (GitHub Container Registry shown):
  ```bash
  echo "$CR_PAT" | docker login ghcr.io -u <user> --password-stdin
  docker push ghcr.io/<org>/cellannot-gpt:<sha>
  ```
- Provision staging infrastructure (e.g., a single VM or ECS task) and run the container:
  ```bash
  docker run -d --name cellannot-staging \
    -e OPENAI_API_KEY=$OPENAI_API_KEY \
    -e ENVIRONMENT=staging \
    -v /srv/cellannot/data:/data \
    -p 8080:8000 \
    ghcr.io/<org>/cellannot-gpt:<sha>
  ```
- Point staging Streamlit/clients to `http://host:8080`.

## 5. Production Deployment
- Mirror staging steps with production secrets and infrastructure (Kubernetes, ECS, or VM).
- Recommended additions:
  - Configure HTTPS termination with reverse proxy (NGINX, Traefik) or managed load balancer.
  - Enable autoscaling with container orchestration (Kubernetes HPA) if throughput spikes are expected.
  - Configure centralized logging/monitoring (e.g., CloudWatch, Grafana Loki, Datadog).

## 6. Secrets Management
- **Local development**: use `.env` (ignored by git) or Streamlit secrets (`.streamlit/secrets.toml`) for UI-specific keys.
- **GitHub Actions**: store `OPENAI_API_KEY`, registry credentials, etc. as Repository Secrets. Reference them via `${{ secrets.OPENAI_API_KEY }}` inside workflows if needed.
- **Runtime (staging/prod)**: prefer secret managers (AWS Secrets Manager, GCP Secret Manager, Hashicorp Vault) or orchestrator secret stores. Mount or export as env vars when starting containers.

## 7. CI/CD Flow
1. GitHub Actions workflow `.github/workflows/ci.yml` runs lint, tests, and Docker build on every commit/PR.
2. Extend workflow with deployment jobs conditioned on branch tags (`main` → production deploy, `develop` → staging deploy).
3. Optionally add image vulnerability scanning (e.g., Trivy) before publishing.

## 8. Troubleshooting
- **Docker build fails fetching dependencies**: ensure network access from build host and the Poetry version specified matches the one used locally.
- **API returns 503 on startup**: run `poetry run python scripts/build_marker_db.py` locally to populate `data/processed` before starting containers, or bake the DB into an object storage bucket the container can fetch at boot.
- **429 / Rate limit errors**: set `OPENAI_REQUESTS_PER_MINUTE` env var to a safe value in production and configure Redis caching if necessary.
