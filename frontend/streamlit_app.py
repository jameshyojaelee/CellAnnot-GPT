from __future__ import annotations

import json
from io import BytesIO
from typing import Any, Dict, List

import altair as alt
import pandas as pd
import streamlit as st
try:
    from fpdf import FPDF
except ImportError:  # pragma: no cover - optional dependency
    FPDF = None  # type: ignore

from frontend.utils import annotate_batch, annotate_cluster_api, format_summary, get_health, status_badge

st.set_page_config(page_title="CellAnnot-GPT", layout="wide")


def _ensure_session_state() -> None:
    st.session_state.setdefault("uploaded_clusters", [])
    st.session_state.setdefault("batch_result", None)
    st.session_state.setdefault("single_cluster_result", None)


def _show_toast(message: str, icon: str = "✅") -> None:
    if hasattr(st, "toast"):
        st.toast(message, icon=icon)


def _render_header(current_page: str, pages: List[str]) -> None:
    st.markdown(
        """
        <div style="display:flex;align-items:center;justify-content:space-between;padding:0.5rem 0;">
          <div style="font-size:1.4rem;font-weight:700;color:#1F6FEB;">CellAnnot-GPT Dashboard</div>
          <div style="font-size:0.9rem;color:#4F6170;">Know your cells in a single click</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    index = pages.index(current_page)
    breadcrumb = " › ".join(pages[: index + 1])
    st.markdown(f"**{breadcrumb}**")
    if len(pages) > 1:
        st.progress(index / (len(pages) - 1))


def _parse_marker_field(value: Any) -> List[str]:
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


def _build_pdf_report(report: Dict[str, Any]) -> bytes:
    if FPDF is None:
        return b"PDF generation requires the fpdf2 package."

    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(0, 10, "CellAnnot-GPT Batch Report", ln=1)

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
        pdf.multi_cell(
            0,
            5,
            f"#{cluster['cluster_id']} | {cluster.get('annotation', {}).get('primary_label', 'Unknown')}"
            f" | status: {cluster.get('status', 'n/a')}\nWarnings: {', '.join(cluster.get('warnings', [])) or 'None'}\n",
        )

    buffer = BytesIO()
    pdf.output(buffer)
    buffer.seek(0)
    return buffer.read()


def load_markers_from_csv(file) -> List[Dict[str, Any]]:
    df = pd.read_csv(file)
    required_cols = {"cluster_id", "markers"}
    if missing := required_cols - set(df.columns):
        raise ValueError(f"CSV missing columns: {', '.join(sorted(missing))}")

    clusters: List[Dict[str, Any]] = []
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
        except Exception as exc:  # noqa: BLE001 - show friendly error
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
            except Exception as exc:  # noqa: BLE001
                st.error(f"Batch annotation failed: {exc}")
                return
        st.session_state["batch_result"] = result
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
    col1.caption("ℹ️ Share of clusters validated successfully.")
    col2.metric(
        "Flagged",
        summary.get("flagged_clusters", 0),
        f"{metrics.get('flagged_rate', 0) * 100:.1f}%",
    )
    col2.caption("ℹ️ Requires manual review or marker cross-checking.")
    col3.metric(
        "Unknown",
        len(summary.get("unknown_clusters", [])),
        f"{metrics.get('unknown_rate', 0) * 100:.1f}%",
    )
    col3.caption("ℹ️ Clusters where the model suggested novel/uncertain identity.")

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
            [{"Reason": key.replace("_", " ").title(), "Count": val} for key, val in flagged_reasons.items()]
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
        file_name="cellannot_report.csv",
        mime="text/csv",
    )
    dl_col2.download_button(
        "⬇️ Download JSON",
        json_data,
        file_name="cellannot_report.json",
        mime="application/json",
    )
    dl_col3.download_button(
        "⬇️ Download PDF",
        pdf_data,
        file_name="cellannot_report.pdf",
        mime="application/pdf",
    )

    st.info(
        "Tip: use your browser's *Print → Save as PDF* or screenshot utilities to capture charts for slides."  # noqa: E501
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

        with st.expander(f"Cluster {cluster_id} → {primary_label} ({status})", expanded=status != "Supported"):
            if markers:
                st.caption(f"Markers: {', '.join(markers)}")
            if confidence:
                st.markdown(f"**Confidence:** {confidence}")
            if warnings:
                st.warning("\n".join(warnings))
            else:
                st.success("Supported by marker database")

            cols = st.columns(2)
            with cols[0]:
                st.write("Annotation")
                st.json(annotation)
            with cols[1]:
                st.write("Validation")
                st.json(validation or {})


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
            payload = {
                "cluster": {"cluster_id": "single", "markers": markers},
                "dataset_context": {k: v for k, v in [("species", species), ("tissue", tissue)] if v},
                "return_validated": True,
            }
            with st.spinner("Annotating..."):
                try:
                    response = annotate_cluster_api(payload)
                except Exception as exc:  # noqa: BLE001
                    st.error(f"Annotation failed: {exc}")
                    return
            st.session_state["single_cluster_result"] = response["result"]
            _show_toast("Single cluster annotated", icon="🔬")

    result = st.session_state.get("single_cluster_result")
    if not result:
        return

    summary = result.get("summary", {})
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
        st.info("Follow-up: inspect markers manually, consult domain experts, or gather additional data.")

    if warnings:
        st.error("\n".join(warnings))
    else:
        st.success("Markers align with known references.")

    st.write("### Annotation Details")
    st.json(annotation)

    st.write("### Validation")
    st.json(validation or {})

    explanation = annotation.get("rationale", "")
    btn = st.download_button(
        "Download Explanation",
        explanation.encode("utf-8"),
        file_name="cluster_explanation.txt",
        mime="text/plain",
    )

def main() -> None:
    _ensure_session_state()

    pages = ["Upload Markers", "Batch Annotate", "Review Results", "Single Cluster"]
    st.sidebar.title("CellAnnot-GPT")
    llm_mode = "unknown"
    badge_label = "API"
    cache_enabled = False
    try:
        health = get_health()
        status = health.get("status", "ok")
        llm_mode = health.get("llm_mode", "unknown")
        cache_enabled = bool(health.get("cache_enabled", False))
        badge_color = "ok" if status == "ok" else "warn"
    except Exception:  # noqa: BLE001 - keep sidebar resilient
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
            "Mock annotator active: predictions use marker heuristics only. Expect reduced accuracy "
            "until an OpenAI API key is provided."
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
    else:
        page_single_cluster()


if __name__ == "__main__":
    main()
