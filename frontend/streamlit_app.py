from __future__ import annotations

import json
from datetime import datetime
from io import BytesIO
from typing import Any

import altair as alt
import pandas as pd
import streamlit as st

try:
    from fpdf import FPDF
except ImportError:  # pragma: no cover - optional dependency
    FPDF = None  # type: ignore

from frontend.utils import (
    annotate_batch,
    annotate_cluster_api,
    build_call_to_action,
    collect_marker_highlights,
    format_marker_links,
    format_summary,
    get_health,
    status_badge,
)

st.set_page_config(page_title="GPT Cell Annotator", layout="wide")


def _ensure_session_state() -> None:
    st.session_state.setdefault("uploaded_clusters", [])
    st.session_state.setdefault("batch_result", None)
    st.session_state.setdefault("single_cluster_result", None)
    st.session_state.setdefault("batch_history", [])


def _show_toast(message: str, icon: str = "✅") -> None:
    if hasattr(st, "toast"):
        st.toast(message, icon=icon)


def _render_header(current_page: str, pages: list[str]) -> None:
    st.markdown(
        """
        <div style="display:flex;align-items:center;justify-content:space-between;
        padding:0.5rem 0;">
          <div style="font-size:1.4rem;font-weight:700;color:#1F6FEB;">
            GPT Cell Annotator Dashboard
          </div>
          <div style="font-size:0.9rem;color:#4F6170;">
            Know your cells in a single click
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    index = pages.index(current_page)
    breadcrumb = " > ".join(pages[: index + 1])
    st.markdown(f"**{breadcrumb}**")
    if len(pages) > 1:
        st.progress(index / (len(pages) - 1))


def _parse_marker_field(value: Any) -> list[str]:
    if isinstance(value, list):
        return [str(v).strip() for v in value if str(v).strip()]
    if isinstance(value, str):
        candidate = value.strip()
        if not candidate:
            return []
        try:
            loaded = json.loads(candidate)
            if isinstance(loaded, list):
                return [str(v).strip() for v in loaded if str(v).strip()]
        except json.JSONDecodeError:
            pass
        # fallback: comma separated string
        return [gene.strip() for gene in candidate.split(",") if gene.strip()]
    return []


def _build_pdf_report(report: dict[str, Any]) -> bytes:
    if FPDF is None:
        return b"PDF generation requires the fpdf2 package."

    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(0, 10, "GPT Cell Annotator Batch Report", ln=1)

    summary = report.get("summary", {})
    metrics = report.get("metrics", {})

    pdf.set_font("Helvetica", size=12)
    pdf.multi_cell(
        0,
        8,
        f"Total clusters: {summary.get('total_clusters', 0)}\n"
        f"Supported: {summary.get('supported_clusters', 0)}\n"
        f"Flagged: {summary.get('flagged_clusters', 0)}\n"
        f"Unknown: {len(summary.get('unknown_clusters', []))}\n",
    )

    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 10, "Metrics", ln=1)
    pdf.set_font("Helvetica", size=11)
    pdf.multi_cell(
        0,
        6,
        f"Support rate: {metrics.get('support_rate', 0)*100:.1f}%\n"
        f"Flagged rate: {metrics.get('flagged_rate', 0)*100:.1f}%\n"
        f"Unknown rate: {metrics.get('unknown_rate', 0)*100:.1f}%",
    )

    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 10, "Clusters", ln=1)
    pdf.set_font("Helvetica", size=10)
    for cluster in report.get("clusters", []):
        primary = cluster.get("annotation", {}).get("primary_label", "Unknown")
        warnings = ", ".join(cluster.get("warnings", [])) or "None"
        content = (
            f"#{cluster['cluster_id']} | {primary} | status: {cluster.get('status', 'n/a')}\n"
            f"Warnings: {warnings}\n"
        )
        pdf.multi_cell(0, 5, content)

    buffer = BytesIO()
    pdf.output(buffer)
    buffer.seek(0)
    return buffer.read()


def load_markers_from_csv(file) -> list[dict[str, Any]]:
    df = pd.read_csv(file)
    required_cols = {"cluster_id", "markers"}
    if missing := required_cols - set(df.columns):
        raise ValueError(f"CSV missing columns: {', '.join(sorted(missing))}")

    clusters: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        clusters.append(
            {
                "cluster_id": str(row["cluster_id"]),
                "markers": _parse_marker_field(row["markers"]),
            }
        )
    return clusters


def page_upload_markers() -> None:
    st.header("Upload Markers")
    st.write(
        "Upload a CSV with columns `cluster_id` and `markers` (JSON array or comma-separated "
        "string). These will be cached for the session."
    )
    uploaded_file = st.file_uploader("Marker CSV", type=["csv"])

    if uploaded_file:
        try:
            clusters = load_markers_from_csv(uploaded_file)
        except Exception as exc:  # - show friendly error
            st.error(f"Failed to parse markers: {exc}")
            return
        st.success(f"Loaded {len(clusters)} clusters.")
        st.json(clusters[:5])
        st.session_state["uploaded_clusters"] = clusters


def page_batch_annotate() -> None:
    st.header("Batch Annotate")
    clusters = st.session_state.get("uploaded_clusters", [])
    if not clusters:
        st.info("Upload marker data first on the Upload Markers page.")
        return

    st.write("Configure dataset context (optional) before running batch annotation.")
    col1, col2 = st.columns(2)
    with col1:
        species = st.text_input("Species", value="Homo sapiens")
    with col2:
        tissue = st.text_input("Tissue / Compartment", value="")

    if st.button("Run Batch Annotation", type="primary"):
        payload = {
            "clusters": clusters,
            "dataset_context": {k: v for k, v in [("species", species), ("tissue", tissue)] if v},
        }
        with st.spinner("Annotating clusters..."):
            try:
                result = annotate_batch(payload)
            except Exception as exc:
                st.error(f"Batch annotation failed: {exc}")
                return
        st.session_state["batch_result"] = result
        history = st.session_state.get("batch_history", [])
        history.append(
            {
                "run_id": datetime.utcnow().strftime("%Y%m%d-%H%M%S"),
                "timestamp": datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),
                "payload": payload,
                "report": result["result"],
            }
        )
        st.session_state["batch_history"] = history[-6:]
        st.success("Batch annotation complete. Review results on the next page.")
        _show_toast("Batch annotation finished", icon="✅")


def page_review_results() -> None:
    st.header("Review Results")
    result = st.session_state.get("batch_result")
    if not result:
        st.info("Run batch annotation to populate this view.")
        return

    report = result["result"]
    summary_str = format_summary(report)
    st.markdown(f"**Summary:** {summary_str}")

    summary = report.get("summary", {})
    metrics = report.get("metrics", {})
    clusters = report.get("clusters", [])

    col1, col2, col3 = st.columns(3)
    col1.metric(
        "Supported",
        summary.get("supported_clusters", 0),
        f"{metrics.get('support_rate', 0) * 100:.1f}%",
    )
    col1.caption("Info: Share of clusters validated successfully.")
    col2.metric(
        "Flagged",
        summary.get("flagged_clusters", 0),
        f"{metrics.get('flagged_rate', 0) * 100:.1f}%",
    )
    col2.caption("Info: Requires manual review or marker cross-checking.")
    col3.metric(
        "Unknown",
        len(summary.get("unknown_clusters", [])),
        f"{metrics.get('unknown_rate', 0) * 100:.1f}%",
    )
    col3.caption("Info: Clusters where the model suggested novel or uncertain identity.")

    chart_data = pd.DataFrame(
        [
            {"Status": "Supported", "Count": summary.get("supported_clusters", 0)},
            {"Status": "Flagged", "Count": summary.get("flagged_clusters", 0)},
        ]
    )
    bar_chart = (
        alt.Chart(chart_data)
        .mark_bar(radius=4)
        .encode(x=alt.X("Status", sort=None), y="Count", color="Status")
    )
    st.altair_chart(bar_chart.interactive(), use_container_width=True)

    confidence_counts = metrics.get("confidence_counts", {})
    if confidence_counts:
        confidence_df = pd.DataFrame(
            [{"Confidence": level, "Count": count} for level, count in confidence_counts.items()]
        )
        conf_chart = (
            alt.Chart(confidence_df)
            .mark_bar(radius=4)
            .encode(x=alt.X("Confidence", sort=None), y="Count", color="Confidence")
        )
        st.altair_chart(conf_chart.interactive(), use_container_width=True)

    flagged_reasons = metrics.get("flagged_reasons", {})
    if flagged_reasons:
        reason_df = pd.DataFrame(
            [
                {
                    "Reason": key.replace("_", " ").title(),
                    "Count": val,
                }
                for key, val in flagged_reasons.items()
            ]
        )
        reason_chart = (
            alt.Chart(reason_df)
            .mark_bar(radius=4)
            .encode(y=alt.Y("Reason", sort="-x"), x="Count", color="Reason")
        )
        st.altair_chart(reason_chart.interactive(), use_container_width=True)

    records = [
        {
            "cluster_id": cluster["cluster_id"],
            "primary_label": cluster.get("annotation", {}).get("primary_label"),
            "status": cluster.get("status"),
            "confidence": cluster.get("confidence"),
            "warnings": " | ".join(cluster.get("warnings", [])),
        }
        for cluster in clusters
    ]
    csv_data = pd.DataFrame(records).to_csv(index=False).encode("utf-8")
    json_data = json.dumps(report, indent=2).encode("utf-8")
    pdf_data = _build_pdf_report(report)

    dl_col1, dl_col2, dl_col3 = st.columns(3)
    dl_col1.download_button(
        "⬇️ Download CSV",
        csv_data,
        file_name="gpt_cell_annotator_report.csv",
        mime="text/csv",
    )
    dl_col2.download_button(
        "⬇️ Download JSON",
        json_data,
        file_name="gpt_cell_annotator_report.json",
        mime="application/json",
    )
    dl_col3.download_button(
        "⬇️ Download PDF",
        pdf_data,
        file_name="gpt_cell_annotator_report.pdf",
        mime="application/pdf",
    )

    st.info(
        "Tip: use your browser's *Print → Save as PDF* or screenshot utilities to "
        "capture charts for slides."
    )

    for cluster in clusters:
        cluster_id = cluster["cluster_id"]
        annotation = cluster.get("annotation", {})
        primary_label = annotation.get("primary_label", "Unknown")
        markers = annotation.get("markers", [])
        warnings = cluster.get("warnings", [])
        validation = cluster.get("validation")
        status = cluster.get("status", "supported").title()
        confidence = cluster.get("confidence")
        highlights = collect_marker_highlights(cluster)
        call_to_action = build_call_to_action(cluster)

        with st.expander(
            f"Cluster {cluster_id} → {primary_label} ({status})",
            expanded=status != "Supported",
        ):
            if markers:
                st.markdown(
                    "**Markers:** " + format_marker_links(markers),
                    unsafe_allow_html=True,
                )
            if confidence:
                st.markdown(f"**Confidence:** {confidence}")
            if warnings:
                st.warning("\n".join(warnings))
            else:
                st.success("Supported by marker database")

            if highlights["warning_markers"]:
                warning_links = format_marker_links(highlights["warning_markers"])
                st.markdown(
                    f"**Markers triggering warnings:** {warning_links}",
                    unsafe_allow_html=True,
                )
            if highlights["supporting_markers"]:
                support_links = format_marker_links(highlights["supporting_markers"])
                st.markdown(
                    f"**Knowledge base support:** {support_links}",
                    unsafe_allow_html=True,
                )

            if call_to_action["next_experiments"] or call_to_action["markers_to_validate"]:
                st.subheader("Next Steps")
                cta_col1, cta_col2 = st.columns(2)
                with cta_col1:
                    st.markdown("**Next experiments**")
                    for item in call_to_action["next_experiments"]:
                        st.write(f"- {item}")
                with cta_col2:
                    st.markdown("**Markers to validate**")
                    if call_to_action["markers_to_validate"]:
                        st.markdown(
                            format_marker_links(call_to_action["markers_to_validate"]),
                            unsafe_allow_html=True,
                        )
                    else:
                        st.write("- None")

            cols = st.columns(2)
            with cols[0]:
                st.write("Annotation")
                st.json(annotation)
            with cols[1]:
                st.write("Validation")
                st.json(validation or {})


def _prepare_run_dataframe(report: dict[str, Any], run_label: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cluster in report.get("clusters", []):
        validation = cluster.get("validation") or {}
        contradictory = validation.get("contradictory_markers") or {}
        rows.append(
            {
                "run": run_label,
                "cluster_id": str(cluster.get("cluster_id")),
                "primary_label": cluster.get("annotation", {}).get("primary_label"),
                "status": cluster.get("status"),
                "contradictions": len(contradictory),
                "warnings": len(cluster.get("warnings") or []),
            }
        )
    return pd.DataFrame(rows)


def page_compare_runs() -> None:
    st.header("Compare Batch Runs")
    history = st.session_state.get("batch_history", [])
    if len(history) < 2:
        st.info("Run at least two batch annotations to unlock comparison mode.")
        return

    options = {run["run_id"]: run for run in history}
    latest_id = history[-1]["run_id"]

    run_a_id = st.selectbox(
        "Primary run",
        list(options.keys()),
        index=list(options.keys()).index(latest_id),
    )
    baseline_choices = [rid for rid in options if rid != run_a_id]
    run_b_id = st.selectbox("Baseline run", baseline_choices, index=0)

    run_a = options[run_a_id]
    run_b = options[run_b_id]

    st.caption(
        f"Comparing **{run_a_id}** ({run_a['timestamp']}) vs **{run_b_id}** ({run_b['timestamp']})"
    )

    df_a = _prepare_run_dataframe(run_a["report"], run_a_id)
    df_b = _prepare_run_dataframe(run_b["report"], run_b_id)
    combined = pd.concat([df_a, df_b], ignore_index=True)

    contradiction_chart = (
        alt.Chart(combined)
        .mark_bar()
        .encode(
            x=alt.X("cluster_id:N", title="Cluster ID", sort=None),
            y=alt.Y("contradictions:Q", title="# Contradictory markers"),
            color=alt.Color("run:N", title="Run"),
            tooltip=["run", "cluster_id", "status", "contradictions", "warnings"],
        )
        .properties(height=320)
    )
    st.altair_chart(contradiction_chart, use_container_width=True)

    pivot = pd.merge(
        df_a,
        df_b,
        on="cluster_id",
        how="outer",
        suffixes=("_new", "_baseline"),
    )
    pivot = pivot.fillna({"status_new": "n/a", "status_baseline": "n/a"})
    columns_to_normalise = [
        "contradictions_new",
        "contradictions_baseline",
        "warnings_new",
        "warnings_baseline",
    ]
    for col in columns_to_normalise:
        if col in pivot.columns:
            pivot[col] = pivot[col].fillna(0).astype(int)
    pivot["status_changed"] = pivot["status_new"] != pivot["status_baseline"]

    columns_in_order = [
        "cluster_id",
        "status_baseline",
        "status_new",
        "primary_label_baseline",
        "primary_label_new",
        "warnings_baseline",
        "warnings_new",
        "contradictions_baseline",
        "contradictions_new",
    ]
    available_columns = [col for col in columns_in_order if col in pivot.columns]
    display_df = pivot[available_columns].rename(
        columns={
            "status_baseline": f"Status {run_b_id}",
            "status_new": f"Status {run_a_id}",
            "primary_label_baseline": f"Label {run_b_id}",
            "primary_label_new": f"Label {run_a_id}",
            "warnings_baseline": f"Warnings {run_b_id}",
            "warnings_new": f"Warnings {run_a_id}",
            "contradictions_baseline": f"Contradictions {run_b_id}",
            "contradictions_new": f"Contradictions {run_a_id}",
        }
    )

    st.subheader("Status Changes")
    st.dataframe(display_df, use_container_width=True)

    highlight_clusters = pivot[
        pivot["status_changed"] | (pivot["contradictions_new"] != pivot["contradictions_baseline"])
    ]
    if not highlight_clusters.empty:
        st.subheader("Clusters to Review")
        for _, row in highlight_clusters.iterrows():
            st.markdown(
                f"- **Cluster {row['cluster_id']}**: contradictions changed from "
                f"{int(row['contradictions_baseline'])} → {int(row['contradictions_new'])}, "
                f"status {row[f'Status {run_b_id}']} → {row[f'Status {run_a_id}']}"
            )


def page_single_cluster() -> None:
    st.header("Single Cluster Annotation")
    species = st.text_input("Species", value="Homo sapiens")
    tissue = st.text_input("Tissue / Compartment", value="")
    markers_text = st.text_area(
        "Marker genes",
        placeholder="Enter markers separated by commas (e.g. MS4A1, CD79A, CD74)",
        height=100,
    )

    markers = [gene.strip() for gene in markers_text.split(",") if gene.strip()]

    if st.button("Annotate Cluster", type="primary"):
        if not markers:
            st.warning("Please provide at least one marker gene.")
        else:
            dataset_context = {}
            if species:
                dataset_context["species"] = species
            if tissue:
                dataset_context["tissue"] = tissue
            payload = {
                "cluster": {"cluster_id": "single", "markers": markers},
                "dataset_context": dataset_context,
                "return_validated": True,
            }
            with st.spinner("Annotating..."):
                try:
                    response = annotate_cluster_api(payload)
                except Exception as exc:
                    st.error(f"Annotation failed: {exc}")
                    return
            st.session_state["single_cluster_result"] = response["result"]
            _show_toast("Single cluster annotated", icon="🔬")

    result = st.session_state.get("single_cluster_result")
    if not result:
        return

    cluster = result.get("clusters", [{}])[0]
    annotation = cluster.get("annotation", {})
    validation = cluster.get("validation")
    warnings = cluster.get("warnings", [])
    status = cluster.get("status", "supported").title()
    confidence = annotation.get("confidence", "Unknown")
    primary_label = annotation.get("primary_label", "Unknown or Novel")

    badge_color = "green" if status == "Supported" else "orange"
    st.markdown(
        f"<span style='padding:6px 12px;border-radius:12px;background:{badge_color};color:white;'>"
        f"{status}</span>",
        unsafe_allow_html=True,
    )
    st.subheader(f"Prediction: {primary_label}")
    st.write(f"**Confidence:** {confidence}")

    if primary_label.lower().startswith("unknown"):
        st.warning("Model suggests this cluster may represent a novel or uncertain cell type.")
        st.info(
            "Follow-up: inspect markers manually, consult domain experts, or gather "
            "additional data."
        )

    if warnings:
        st.error("\n".join(warnings))
    else:
        st.success("Markers align with known references.")

    st.write("### Annotation Details")
    st.json(annotation)

    st.write("### Validation")
    st.json(validation or {})

    explanation = annotation.get("rationale", "")
    st.download_button(
        "Download Explanation",
        explanation.encode("utf-8"),
        file_name="cluster_explanation.txt",
        mime="text/plain",
    )


def main() -> None:
    _ensure_session_state()

    pages = [
        "Upload Markers",
        "Batch Annotate",
        "Review Results",
        "Compare Batches",
        "Single Cluster",
    ]
    st.sidebar.title("GPT Cell Annotator")
    st.sidebar.markdown(
        "[Getting Started](docs/getting_started.md) · "
        "[API Reference](docs/api_reference.md) · "
        "[Operations](docs/operations.md)",
        unsafe_allow_html=True,
    )
    llm_mode = "unknown"
    badge_label = "API"
    cache_enabled = False
    try:
        health = get_health()
        status = health.get("status", "ok")
        llm_mode = health.get("llm_mode", "unknown")
        cache_enabled = bool(health.get("cache_enabled", False))
        badge_color = "ok" if status == "ok" else "warn"
    except Exception:  # keep sidebar resilient
        badge_color = "error"
    else:
        if llm_mode == "mock":
            badge_label = "API (mock)"
            badge_color = "warn"
    st.sidebar.markdown(status_badge(badge_label, badge_color), unsafe_allow_html=True)
    if llm_mode == "mock":
        st.sidebar.warning(
            "LLM requests are running in mock mode. Outputs use heuristic rules; configure "
            "`OPENAI_API_KEY` for live annotations."
        )
        st.info(
            "Mock annotator active: predictions use marker heuristics only."
            " Expect reduced accuracy until an OpenAI API key is provided."
        )
    if cache_enabled:
        st.sidebar.success("Redis cache active for repeated annotations.")

    page = st.sidebar.radio("Navigation", pages)

    _render_header(page, pages)

    if page == "Upload Markers":
        page_upload_markers()
    elif page == "Batch Annotate":
        page_batch_annotate()
    elif page == "Review Results":
        page_review_results()
    elif page == "Compare Batches":
        page_compare_runs()
    else:
        page_single_cluster()


if __name__ == "__main__":
    main()
