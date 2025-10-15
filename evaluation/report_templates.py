from __future__ import annotations

from collections.abc import Iterable


def _format_delta(delta: float | None, regression_threshold: float) -> str:
    if delta is None:
        return ""
    delta_percent = f"{delta:+.2%}"
    if delta <= -regression_threshold:
        return f" **{delta_percent}**"
    return f" {delta_percent}"


def render_dataset_report(
    dataset_summary: dict[str, object],
    *,
    deltas: dict[str, dict[str, float | None]] | None = None,
    regression_threshold: float = 0.05,
) -> str:
    lines = [
        f"# Benchmark Report – {dataset_summary.get('dataset', 'unknown')}",
        "",
    ]
    for model in dataset_summary.get("models", []):
        name = model["model_name"]
        model_deltas = (deltas or {}).get(name, {})
        lines.extend(
            [
                f"## {name}",
                "",
                (
                    f"- Accuracy: {model['accuracy']:.2%}"
                    f"{_format_delta(model_deltas.get('accuracy'), regression_threshold)}"
                ),
                (
                    f"- Macro F1: {model['macro_f1']:.2%}"
                    f"{_format_delta(model_deltas.get('macro_f1'), regression_threshold)}"
                ),
                f"- Top-3 recall: {model['top3_recall']:.2%}",
                f"- Unknown rate: {model['unknown_rate']:.2%}",
                f"- Flag precision: {model['flag_precision']:.2%}",
                f"- Time per cluster: {model['time_per_cluster']:.3f}s",
                "",
                "### Per-class metrics",
            ]
        )
        for label, metrics in model["per_class"].items():
            lines.append(
                "- **{label}**: Precision {precision:.2%}, Recall {recall:.2%}, "
                "F1 {f1:.2%}".format(
                    label=label,
                    precision=metrics["precision"],
                    recall=metrics["recall"],
                    f1=metrics["f1"],
                )
            )
        lines.append("")
        lines.append("### Prediction details")
        for pred in model["predictions"]:
            lines.append(
                "- Cluster {id}: predicted {pred} (truth {truth}); status {status}".format(
                    id=pred["cluster_id"],
                    pred=pred["primary_label"],
                    truth=pred["ground_truth"],
                    status=pred.get("status", "unknown"),
                )
            )
        lines.append("")
    return "\n".join(lines)


def render_sparkline_csv(dataset: str, history: Iterable[dict[str, float]]) -> str:
    rows: list[str] = ["dataset,date,accuracy,macro_f1"]
    for entry in history:
        rows.append(
            f"{dataset},{entry.get('date', '')},"
            f"{entry.get('accuracy', 0.0):.4f},"
            f"{entry.get('macro_f1', 0.0):.4f}"
        )
    return "\n".join(rows)


__all__ = ["render_dataset_report", "render_sparkline_csv"]
