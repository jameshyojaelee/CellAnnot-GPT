from __future__ import annotations

import json
from typing import Any, Dict, List

import altair as alt
import pandas as pd
import streamlit as st

from frontend.utils import annotate_batch, format_summary, get_health, status_badge

st.set_page_config(page_title="CellAnnot-GPT", layout="wide")


def _ensure_session_state() -> None:
    st.session_state.setdefault("uploaded_clusters", [])
    st.session_state.setdefault("batch_result", None)


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

    if st.button("Run Batch Annotation"):
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
    col2.metric(
        "Flagged",
        summary.get("flagged_clusters", 0),
        f"{metrics.get('flagged_rate', 0) * 100:.1f}%",
    )
    col3.metric(
        "Unknown",
        len(summary.get("unknown_clusters", [])),
        f"{metrics.get('unknown_rate', 0) * 100:.1f}%",
    )

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
    st.altair_chart(bar_chart, use_container_width=True)

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
        st.altair_chart(conf_chart, use_container_width=True)

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
        st.altair_chart(reason_chart, use_container_width=True)

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


def main() -> None:
    _ensure_session_state()

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

    page = st.sidebar.radio(
        "Navigation",
        ("Upload Markers", "Batch Annotate", "Review Results"),
    )

    if page == "Upload Markers":
        page_upload_markers()
    elif page == "Batch Annotate":
        page_batch_annotate()
    else:
        page_review_results()


if __name__ == "__main__":
    main()
