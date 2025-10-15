from __future__ import annotations

from collections.abc import Iterable


def _format_delta(delta: float | None, regression_threshold: float) -> str:
    if delta is None:
        return ""
    delta_percent = f"{delta:+.2%}"
    if delta <= -regression_threshold:
        return f" **{delta_percent}**"
    return f" {delta_percent}"


def render_markdown_report(
    result: dict,
    *,
    deltas: dict[str, float | None] | None = None,
    regression_threshold: float = 0.05,
) -> str:
    lines = [
        "# CellAnnot-GPT Benchmark Report",
        "",
        (
            f"- Accuracy: {result['accuracy']:.2%}"
            f"{_format_delta((deltas or {}).get('accuracy'), regression_threshold)}"
        ),
        (
            f"- Macro F1: {result['macro_f1']:.2%}"
            f"{_format_delta((deltas or {}).get('macro_f1'), regression_threshold)}"
        ),
        "",
        "## Per-class Metrics",
        "",
    ]
    for label, metrics in result["per_class"].items():
        lines.append(
            "- **{label}**: Precision {precision:.2%}, "
            "Recall {recall:.2%}, F1 {f1:.2%}".format(
                label=label,
                precision=metrics["precision"],
                recall=metrics["recall"],
                f1=metrics["f1"],
            )
        )
    lines.append("\n## Prediction Details\n")
    for pred in result["predictions"]:
        lines.append(
            "- Cluster {id}: predicted {pred} (truth {truth})".format(
                id=pred["cluster_id"],
                pred=pred["predicted"],
                truth=pred["ground_truth"],
            )
        )
    return "\n".join(lines)


def render_text_confusion_matrix(result: dict) -> str:
    truth_labels = sorted({row["ground_truth"] for row in result["predictions"]})
    pred_labels = sorted({row["predicted"] for row in result["predictions"]})
    labels = sorted(set(truth_labels) | set(pred_labels))
    matrix_lines = ["Confusion Matrix (rows=truth, cols=prediction)"]
    header = "\t".join([" ", *labels])
    matrix_lines.append(header)
    counts = {
        label: {label_name: 0 for label_name in labels}
        for label in labels
    }
    for row in result["predictions"]:
        truth = row.get("ground_truth")
        pred = row.get("predicted")
        if truth not in counts:
            counts[truth] = {label_name: 0 for label_name in labels}
        if pred not in counts:
            for key in counts:
                counts[key].setdefault(pred, 0)
        counts[truth][pred] += 1
    for truth in labels:
        row_counts = "\t".join(str(counts[truth][pred]) for pred in labels)
        matrix_lines.append(f"{truth}\t{row_counts}")
    return "\n".join(matrix_lines)


def render_sparkline_csv(dataset: str, history: Iterable[dict[str, float]]) -> str:
    rows: list[str] = ["dataset,date,accuracy,macro_f1"]
    for entry in history:
        rows.append(
            f"{dataset},{entry.get('date', '')},"
            f"{entry.get('accuracy', 0.0):.4f},"
            f"{entry.get('macro_f1', 0.0):.4f}"
        )
    return "\n".join(rows)
