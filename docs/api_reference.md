# API Reference

GPT Cell Annotator exposes a FastAPI backend. Below is a summary of the main endpoints.

## `GET /health`
- Returns service status, LLM mode (`live` or `mock`), and cache status.

Example response:
```json
{
  "status": "ok",
  "llm_mode": "live",
  "cache_enabled": true
}
```

## `POST /annotate_cluster`
Annotate a single cluster. Payload supports an optional `return_validated` flag to include validation report.

Request body:
```json
{
  "cluster": {
    "cluster_id": "0",
    "markers": ["MS4A1", "CD79A"]
  },
  "dataset_context": {
    "species": "Homo sapiens",
    "tissue": "Blood"
  },
  "return_validated": true
}
```

Response body:
```json
{
  "result": {
    "summary": {
      "total_clusters": 1,
      "supported_clusters": 1,
      "flagged_clusters": 0,
      "unknown_clusters": []
    },
    "metrics": {
      "support_rate": 1.0,
      "flagged_rate": 0.0,
      "unknown_rate": 0.0,
      "flagged_reasons": {},
      "confidence_counts": {"High": 1}
    },
    "clusters": [
      {
        "cluster_id": "0",
        "status": "supported",
        "confidence": "High",
        "warnings": [],
        "annotation": {
          "primary_label": "B cell",
          "ontology_id": "CL:0000236",
          "confidence": "High",
          "rationale": "Matched markers: MS4A1, CD79A"
        },
        "validation": {
          "cluster_id": "0",
          "primary_label": "B cell",
          "is_supported": true
        }
      }
    ]
  }
}
```

## `POST /annotate_batch`
Annotate multiple clusters in one request.

Request body:
```json
{
  "clusters": [
    {"cluster_id": "0", "markers": ["MS4A1"]},
    {"cluster_id": "1", "markers": ["CD3E"]}
  ],
  "dataset_context": {"species": "Homo sapiens"}
}
```

Response body includes dataset summary, metrics, and per-cluster breakdown analogous to the single-cluster endpoint.

## OpenAPI Schema
Run the backend locally and visit `http://localhost:8000/docs` for interactive Swagger UI or download OpenAPI JSON from `http://localhost:8000/openapi.json`.

## Authentication
The current prototype does not require auth. For production, consider adding API keys or OAuth headers.
