"""Streamlit UI for session-based exploratory analysis.

Launch with: ``streamlit run ui/app.py``
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import streamlit as st
import polars as pl

from ui.session import Session
from ui.registry import AssetRegistry
from ui import plotting

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

CONFIG_PATH = os.environ.get("UI_CONFIG_PATH", "ui/config.json")
ALL_STATS = ["enrichment", "rate_ratio", "auc", "auprc", "all"]


def _get_session() -> Session:
    """Get or create the Streamlit-cached Session."""
    if "session" not in st.session_state:
        name = st.session_state.get("session_name", "default")
        st.session_state["session"] = Session(
            config_path=CONFIG_PATH,
            session_name=name,
        )
    return st.session_state["session"]


def _switch_session(name: str) -> None:
    """Close current session and open a new one."""
    if "session" in st.session_state:
        st.session_state["session"].close()
    st.session_state["session_name"] = name
    st.session_state["session"] = Session(
        config_path=CONFIG_PATH,
        session_name=name,
    )


def main() -> None:
    st.set_page_config(
        page_title="Genetics Gym Explorer",
        page_icon="🧬",
        layout="wide",
    )
    st.title("Genetics Gym Explorer")

    session = _get_session()
    registry = session.registry

    # ------------------------------------------------------------------ #
    # Sidebar: Session + Asset Browser
    # ------------------------------------------------------------------ #
    with st.sidebar:
        _render_session_picker(session)
        st.divider()
        _render_score_selector(session, registry)
        st.divider()
        _render_eval_selector(session, registry)
        st.divider()
        _render_filter_selector(session, registry)
        st.divider()
        _render_pipeline_settings(session, registry)

    # ------------------------------------------------------------------ #
    # Main Panel
    # ------------------------------------------------------------------ #
    col_left, col_right = st.columns([1, 2])

    with col_left:
        _render_active_scores(session)
        _render_analysis_controls(session)

    with col_right:
        _render_results()

    # ------------------------------------------------------------------ #
    # Session Log (collapsible)
    # ------------------------------------------------------------------ #
    with st.expander("Session Log", expanded=False):
        _render_session_log(session)


# ------------------------------------------------------------------ #
# Sidebar Components
# ------------------------------------------------------------------ #

def _render_session_picker(session: Session) -> None:
    st.subheader("Session")
    existing = session.engine.list_cached_sessions()
    if not existing:
        existing = ["default"]

    current_idx = 0
    current_name = session.session_name
    if current_name in existing:
        current_idx = existing.index(current_name)

    selected = st.selectbox("Active session", existing, index=current_idx, key="sb_session")
    if selected and selected != current_name:
        _switch_session(selected)
        st.rerun()

    col1, col2 = st.columns(2)
    with col1:
        new_name = st.text_input("New session name", key="new_session_name", label_visibility="collapsed", placeholder="new session name")
    with col2:
        if st.button("Create", key="btn_new_session"):
            if new_name and new_name.strip():
                _switch_session(new_name.strip())
                st.rerun()

    if st.button("Rebuild Session", key="btn_rebuild"):
        with st.spinner("Rebuilding..."):
            session.rebuild()
        st.success("Session rebuilt from scratch")


def _render_score_selector(session: Session, registry: AssetRegistry) -> None:
    st.subheader("Prediction Tables")

    variant_preds = registry.variant_predictions()
    gene_preds = registry.gene_predictions()

    if variant_preds:
        st.caption("Variant-level")
        all_variant_scores = []
        for entry in variant_preds:
            all_variant_scores.extend(entry.score_columns)

        selected_scores = st.multiselect(
            "Variant scores",
            options=all_variant_scores,
            default=session.score_cols if session.score_cols else [],
            key="ms_variant_scores",
            label_visibility="collapsed",
        )
        st.session_state["selected_variant_scores"] = selected_scores

    if gene_preds:
        st.caption("Gene-level")
        all_gene_scores = []
        for entry in gene_preds:
            all_gene_scores.extend(entry.score_columns)

        selected_gene = st.multiselect(
            "Gene scores",
            options=all_gene_scores,
            default=[],
            key="ms_gene_scores",
            label_visibility="collapsed",
        )
        st.session_state["selected_gene_scores"] = selected_gene

    if st.button("Load Selected Scores", key="btn_load_scores"):
        all_selected = (
            st.session_state.get("selected_variant_scores", [])
            + st.session_state.get("selected_gene_scores", [])
        )
        if all_selected:
            with st.spinner(f"Loading {len(all_selected)} scores..."):
                join_type = st.session_state.get("join_type", "inner")
                session.load_base(all_selected, join_type=join_type)
            st.success(f"Loaded {len(all_selected)} scores")
        else:
            st.warning("No scores selected")


def _render_eval_selector(session: Session, registry: AssetRegistry) -> None:
    st.subheader("Evaluation Tables")

    evals = registry.list_evaluations()
    eval_names = []
    for e in evals:
        stem = os.path.basename(e.path).split(".")[0]
        short = stem.split("-")[0] if "-" in stem else stem
        eval_names.append(short)

    for i, (e, name) in enumerate(zip(evals, eval_names)):
        is_attached = name in session._eval_aliases
        if st.checkbox(
            f"{name} ({e.level})",
            value=is_attached,
            key=f"cb_eval_{i}",
        ):
            if not is_attached:
                with st.spinner(f"Attaching {name}..."):
                    session.attach_eval(name)
        else:
            if is_attached:
                session.detach_eval(name)


def _render_filter_selector(session: Session, registry: AssetRegistry) -> None:
    st.subheader("Filters")
    filters = registry.list_filters()
    if not filters:
        st.caption("No filters configured")
        return

    for i, path in enumerate(filters):
        from ui.prepare import _derive_filter_name
        fname = _derive_filter_name(path)
        short_name = fname.split("_", 1)[-1] if "_" in fname else fname
        display = f"{short_name} ({os.path.basename(os.path.dirname(path))})"

        enabled = st.checkbox(display, value=False, key=f"cb_filter_{i}")
        if enabled:
            session.toggle_filter(fname, True)


def _render_pipeline_settings(session: Session, registry: AssetRegistry) -> None:
    st.subheader("Pipeline Settings")

    join_type = st.selectbox(
        "Join Type",
        ["inner", "left", "full outer"],
        index=0,
        key="join_type",
    )
    st.session_state["join_type"] = join_type


# ------------------------------------------------------------------ #
# Main Panel Components
# ------------------------------------------------------------------ #

def _render_active_scores(session: Session) -> None:
    st.subheader("Active Scores")

    if not session.score_cols:
        st.info("No scores loaded. Select scores in the sidebar and click 'Load Selected Scores'.")
        return

    for col in session.score_cols:
        c1, c2 = st.columns([3, 1])
        with c1:
            st.text(col)
        with c2:
            if st.button("Remove", key=f"btn_rm_{col}"):
                session.remove_score(col)
                st.rerun()

    info = session.engine.table_info()
    if info.get("row_count"):
        st.metric("Rows", f"{info['row_count']:,}")
        st.metric("Columns", len(info.get("columns", {})))


def _render_analysis_controls(session: Session) -> None:
    st.subheader("Analysis")

    stat = st.selectbox("Statistic", ALL_STATS, key="sel_stat")
    eval_set = st.selectbox(
        "Eval Set",
        options=session._eval_aliases if session._eval_aliases else ["(none attached)"],
        key="sel_eval_set",
    )
    threshold_str = st.text_input(
        "Thresholds",
        value=", ".join(str(t) for t in session.registry.pipeline.percentile_thresholds),
        key="txt_thresholds",
    )

    col1, col2 = st.columns(2)
    with col1:
        if st.button("Run Analysis", key="btn_run", type="primary"):
            if eval_set == "(none attached)":
                st.warning("Attach an eval set first")
                return
            thresholds = [float(t.strip()) for t in threshold_str.split(",") if t.strip()]
            with st.spinner("Running statistics..."):
                try:
                    results, timings, vsm_cmp = session.run_stat(
                        stat=stat,
                        eval_set=eval_set,
                        thresholds=thresholds,
                    )
                    st.session_state["results"] = results
                    st.session_state["timings"] = timings
                    st.session_state["vsm_cmp"] = vsm_cmp
                    st.success(f"Done! {results.height} result rows")
                except Exception as exc:
                    st.error(f"Analysis failed: {exc}")

    with col2:
        if st.button("Compare Scores", key="btn_compare"):
            if eval_set == "(none attached)":
                st.warning("Attach an eval set first")
                return
            thresholds = [float(t.strip()) for t in threshold_str.split(",") if t.strip()]
            with st.spinner("Comparing..."):
                try:
                    results = session.compare(
                        stat=stat,
                        eval_set=eval_set,
                        thresholds=thresholds,
                    )
                    st.session_state["results"] = results
                    st.success(f"Comparison: {results.height} rows")
                except Exception as exc:
                    st.error(f"Comparison failed: {exc}")

    st.divider()

    col_save, col_cli = st.columns(2)
    with col_save:
        if st.button("Save to Disk", key="btn_save"):
            output = session.registry.output
            with st.spinner(f"Exporting to {output}..."):
                session.save(output)
            st.success(f"Saved: {output}")

    with col_cli:
        if st.button("Export CLI Command", key="btn_cli"):
            cmd = session.export_cli_command()
            st.code(cmd, language="bash")


def _render_results() -> None:
    st.subheader("Results")

    results = st.session_state.get("results")
    if results is None or results.is_empty():
        st.info("Run an analysis to see results here.")
        return

    tab_table, tab_bar, tab_line, tab_heatmap, tab_ci = st.tabs([
        "Table", "Bar Chart", "Line Chart", "Heatmap", "CI Plot",
    ])

    with tab_table:
        st.dataframe(results.to_pandas(), use_container_width=True, height=400)

    with tab_bar:
        eval_names = results.get_column("eval_name").unique().to_list() if "eval_name" in results.columns else []
        eval_sel = eval_names[0] if eval_names else None
        fig = plotting.enrichment_bar(results, eval_name=eval_sel)
        st.plotly_chart(fig, use_container_width=True)

    with tab_line:
        fig = plotting.threshold_line(results, eval_name=eval_sel if eval_names else None)
        st.plotly_chart(fig, use_container_width=True)

    with tab_heatmap:
        fig = plotting.score_eval_heatmap(results)
        st.plotly_chart(fig, use_container_width=True)

    with tab_ci:
        fig = plotting.enrichment_dot_with_ci(results, eval_name=eval_sel if eval_names else None)
        st.plotly_chart(fig, use_container_width=True)


# ------------------------------------------------------------------ #
# Session Log
# ------------------------------------------------------------------ #

def _render_session_log(session: Session) -> None:
    history = session.history
    if not history:
        st.caption("No actions recorded yet.")
        return

    for action in reversed(history):
        ts = action.timestamp.strftime("%H:%M:%S")
        params_str = ", ".join(f"{k}={v}" for k, v in action.params.items() if v)
        st.text(f"{ts}  {action.type}: {params_str}")

    col1, col2 = st.columns(2)
    with col1:
        if st.button("Undo Last", key="btn_undo"):
            session.undo()
            st.rerun()
    with col2:
        db_history = session.engine.get_history(limit=20)
        if db_history:
            st.caption("Persisted history (from DuckDB):")
            for h in db_history:
                st.text(f"  [{h['seq']}] {h['action']}: {h['params']}")


if __name__ == "__main__":
    main()
