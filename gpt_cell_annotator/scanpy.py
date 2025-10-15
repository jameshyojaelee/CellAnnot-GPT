"""Scanpy integration helpers and CLI entrypoints for GPT Cell Annotator."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData

from backend.llm.annotator import Annotator
from backend.validation.crosscheck import crosscheck_batch
from backend.validation.report import DatasetReport, build_structured_report

try:  # Optional dependency; only needed when we must compute marker rankings.
    import scanpy as sc
except ImportError:  # pragma: no cover - handled in annotate_anndata
    sc = None


def _ensure_rankings(
    adata: AnnData,
    cluster_key: str,
    *,
    top_n_markers: int,
    method: str = "wilcoxon",
) -> None:
    """Compute rank_genes_groups if missing or for a different cluster key."""

    rankings = adata.uns.get("rank_genes_groups")
    params = (rankings or {}).get("params", {})
    if rankings is not None and params.get("groupby") == cluster_key:
        return

    if sc is None:
        raise ImportError(
            "scanpy is required to compute marker rankings. Install scanpy or "
            "precompute `rank_genes_groups` for the AnnData object."
        )

    sc.tl.rank_genes_groups(adata, groupby=cluster_key, n_genes=top_n_markers, method=method)


def _unique_ordered(markers: Iterable[Any]) -> list[str]:
    """Return unique marker symbols while preserving their first-seen order."""

    seen: set[str] = set()
    ordered: list[str] = []
    for marker in markers:
        if not isinstance(marker, str):
            continue
        normalised = marker.strip()
        if not normalised:
            continue
        upper = normalised.upper()
        if upper in seen:
            continue
        seen.add(upper)
        ordered.append(upper)
    return ordered


def _extract_markers(adata: AnnData, cluster_key: str, top_n_markers: int) -> dict[str, list[str]]:
    """Derive top markers for each cluster from rank_genes_groups results."""

    rankings = adata.uns.get("rank_genes_groups")
    if rankings is None or "names" not in rankings:
        raise ValueError(
            f"AnnData object is missing `rank_genes_groups` results for cluster key '{cluster_key}'. "
            "Run scanpy.tl.rank_genes_groups beforehand or allow annotate_anndata to compute them."
        )

    names = rankings["names"]
    if isinstance(names, np.ndarray) and names.dtype.names:
        groups = names.dtype.names
        cluster_to_markers: dict[str, list[str]] = {}
        for group in groups:
            markers = _unique_ordered(names[group][:top_n_markers])
            cluster_to_markers[str(group)] = markers
        return cluster_to_markers

    if isinstance(names, Mapping):
        return {
            str(cluster): _unique_ordered(values[:top_n_markers])
            for cluster, values in names.items()
        }

    raise ValueError("Unsupported structure for `rank_genes_groups['names']`.")


def _load_marker_db(marker_db: pd.DataFrame | None, marker_db_path: Path | None) -> pd.DataFrame:
    if marker_db is not None:
        return marker_db[MARKER_DB_COLUMNS]
    path = marker_db_path or Path("data/processed/marker_db.parquet")
    df = pd.read_parquet(path)
    return df[MARKER_DB_COLUMNS]


def annotate_anndata(
    adata: AnnData,
    cluster_key: str,
    *,
    species: str,
    tissue: str | None = None,
    top_n_markers: int = 5,
    result_prefix: str = "gptca",
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    compute_rankings: bool = True,
) -> tuple[AnnData, dict[str, Any]]:
    """Annotate clusters in an AnnData object using GPT Cell Annotator."""

    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not found in adata.obs.")

    if compute_rankings:
        _ensure_rankings(
            adata,
            cluster_key,
            top_n_markers=top_n_markers,
        )

    cluster_markers = _extract_markers(adata, cluster_key, top_n_markers)

    annotator_instance = annotator or Annotator()

    dataset_context: dict[str, str] = {"species": species}
    if tissue:
        dataset_context["tissue"] = tissue

    clusters_payload = [
        {"cluster_id": cluster_id, "markers": markers}
        for cluster_id, markers in cluster_markers.items()
    ]

    raw_result = annotator_instance.annotate_batch(clusters_payload, dataset_context)

    marker_df = _load_marker_db(
        marker_db,
        Path(marker_db_path) if marker_db_path is not None else None,
    )

    annotations: list[dict[str, Any]] = []
    for cluster in clusters_payload:
        cluster_id = cluster["cluster_id"]
        cluster_result = dict(raw_result.get(cluster_id) or {})
        cluster_result.setdefault("primary_label", "Unknown or Novel")
        cluster_result.setdefault("confidence", "Unknown")
        cluster_result.setdefault("rationale", "")
        annotations.append(
            {
                **cluster_result,
                "cluster_id": cluster_id,
                "markers": cluster["markers"],
            }
        )

    crosschecked = crosscheck_batch(
        annotations,
        marker_df,
        species=species,
        tissue=tissue,
    )
    report_model: DatasetReport = build_structured_report(annotations, crosschecked)
    report = report_model.model_dump()

    label_col = f"{result_prefix}_label"
    confidence_col = f"{result_prefix}_confidence"
    status_col = f"{result_prefix}_status"
    rationale_col = f"{result_prefix}_rationale"
    ontology_col = f"{result_prefix}_ontology_id"

    clusters_index = _build_cluster_index(adata.obs[cluster_key].astype(str))
    adata.obs[label_col] = None
    adata.obs[confidence_col] = None
    adata.obs[status_col] = None
    adata.obs[rationale_col] = None
    adata.obs[ontology_col] = None

    for cluster_report in report["clusters"]:
        cluster_id = str(cluster_report["cluster_id"])
        obs_indices = clusters_index.get(cluster_id, [])
        annotation = cluster_report.get("annotation", {})
        label = annotation.get("primary_label")
        confidence = cluster_report.get("confidence")
        status = cluster_report.get("status")
        rationale = annotation.get("rationale")
        ontology = annotation.get("ontology_id")
        for obs_idx in obs_indices:
            adata.obs.at[obs_idx, label_col] = label
            adata.obs.at[obs_idx, confidence_col] = confidence
            adata.obs.at[obs_idx, status_col] = status
            adata.obs.at[obs_idx, rationale_col] = rationale
            adata.obs.at[obs_idx, ontology_col] = ontology

    return adata, report


def _build_cluster_index(cluster_series: pd.Series) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for obs_idx, cluster_id in cluster_series.items():
        mapping.setdefault(cluster_id, []).append(obs_idx)
    return mapping


def _report_to_dataframe(report: dict[str, Any]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for cluster in report.get("clusters", []):
        annotation = cluster.get("annotation", {})
        records.append(
            {
                "cluster_id": cluster.get("cluster_id"),
                "primary_label": annotation.get("primary_label"),
                "confidence": cluster.get("confidence"),
                "status": cluster.get("status"),
                "ontology_id": annotation.get("ontology_id"),
                "rationale": annotation.get("rationale"),
                "warnings": "; ".join(cluster.get("warnings") or []),
            }
        )
    return pd.DataFrame.from_records(records)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m gpt_cell_annotator.scanpy",
        description="Scanpy integration tools for GPT Cell Annotator.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate clusters in an AnnData .h5ad file.",
    )
    annotate_parser.add_argument("input", type=Path, help="Path to the input .h5ad file.")
    annotate_parser.add_argument(
        "--cluster-key",
        required=True,
        help="Column in adata.obs that identifies cluster assignments.",
    )
    annotate_parser.add_argument(
        "--species",
        required=True,
        help="Species context to provide to the annotator (e.g. 'Homo sapiens').",
    )
    annotate_parser.add_argument(
        "--tissue",
        help="Optional tissue context for the dataset.",
    )
    annotate_parser.add_argument(
        "--top-n-markers",
        type=int,
        default=5,
        help="Number of top markers per cluster to send for annotation.",
    )
    annotate_parser.add_argument(
        "--prefix",
        default="gptca",
        help="Prefix for columns written to adata.obs.",
    )
    annotate_parser.add_argument(
        "--output",
        type=Path,
        help="Destination .h5ad path. Defaults to in-place overwrite.",
    )
    annotate_parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Optional path to write a per-cluster summary CSV.",
    )
    annotate_parser.add_argument(
        "--summary-json",
        type=Path,
        help="Optional path to write the full annotation report JSON.",
    )
    annotate_parser.add_argument(
        "--marker-db",
        type=Path,
        help="Path to marker_db.parquet. Defaults to data/processed/marker_db.parquet.",
    )
    annotate_parser.add_argument(
        "--skip-recompute-markers",
        action="store_true",
        help="Assume rank_genes_groups already computed; do not call scanpy.tl.rank_genes_groups.",
    )

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    adata = ad.read_h5ad(args.input)
    _, report = annotate_anndata(
        adata,
        args.cluster_key,
        species=args.species,
        tissue=args.tissue,
        top_n_markers=args.top_n_markers,
        result_prefix=args.prefix,
        marker_db_path=args.marker_db,
        compute_rankings=not args.skip_recompute_markers,
    )

    output_path = args.output or args.input
    adata.write(output_path)

    if args.summary_csv:
        df = _report_to_dataframe(report)
        df.to_csv(args.summary_csv, index=False)
    if args.summary_json:
        args.summary_json.write_text(json.dumps(report, indent=2), encoding="utf-8")

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
MARKER_DB_COLUMNS = [
    "source",
    "cell_type",
    "ontology_id",
    "gene_symbol",
    "species",
    "tissue",
    "evidence",
    "reference",
    "evidence_score",
]
