from __future__ import annotations

from typing import Dict


def render_markdown_report(result: Dict) -> str:
    lines = [
        "# CellAnnot-GPT Benchmark Report",
        "",
        f"- Accuracy: {result['accuracy']:.2%}",
        f"- Macro F1: {result['macro_f1']:.2%}",
        "",
        "## Per-class Metrics",
        "",
    ]
    for label, metrics in result["per_class"].items():
        lines.append(f"- **{label}**: Precision {metrics['precision']:.2%}, Recall {metrics['recall']:.2%}, F1 {metrics['f1']:.2%}")
    lines.append("\n## Prediction Details\n")
    for pred in result["predictions"]:
        lines.append(
            f"- Cluster {pred['cluster_id']}: predicted {pred['predicted']} (truth {pred['ground_truth']})"
        )
    return "\n".join(lines)


def render_text_confusion_matrix(result: Dict) -> str:
    labels = sorted(result["per_class"].keys())
    matrix_lines = ["Confusion Matrix (rows=truth, cols=prediction)"]
    header = "\t".join([" "] + labels)
    matrix_lines.append(header)
    counts = {label: {l: 0 for l in labels} for label in labels}
    for row in result["predictions"]:
        counts[row["ground_truth"]][row["predicted"]] += 1
    for truth in labels:
        row_counts = "\t".join(str(counts[truth][pred]) for pred in labels)
        matrix_lines.append(f"{truth}\t{row_counts}")
    return "\n".join(matrix_lines)
