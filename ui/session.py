"""Session class: composable operations built on the DuckDB engine.

This is the primary API for both the Streamlit UI and programmatic/notebook use.
Each Session corresponds to a named ``.duckdb`` file in ``local_cache/``.

Usage::

    from ui.session import Session

    s = Session(config_path="ui/config.json", session_name="my_analysis")
    s.load_base(["AM_score", "esm1b", "cadd_score"], join_type="inner")
    s.attach_eval("dd")
    results, timings, vsm_cmp = s.run_stat("enrichment", eval_set="dd")
    s.save("output.parquet", include_percentiles=True)
"""

from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable

import polars as pl

from ui.engine import DuckDBEngine
from ui.prepare import Preprocessor
from ui.registry import AssetRegistry, PredictionEntry, EvaluationEntry

logger = logging.getLogger(__name__)


@dataclass
class Action:
    type: str
    params: dict[str, Any]
    timestamp: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    undo_sql: str | None = None


class Session:
    """High-level composable session for interactive data exploration.

    Wraps a :class:`DuckDBEngine` and :class:`AssetRegistry`, providing
    named operations (load, swap, attach, run, save) with undo support
    and session history.
    """

    def __init__(
        self,
        config_path: str = "ui/config.json",
        session_name: str = "default",
        *,
        auto_prepare: bool = True,
    ) -> None:
        self.config_path = config_path
        self.session_name = session_name
        self.registry = AssetRegistry(config_path)

        staging_dir = os.path.join(
            os.path.dirname(config_path), "local_cache", "staging"
        )
        cache_dir = os.path.join(os.path.dirname(config_path), "local_cache")
        manifest_path = os.path.join(staging_dir, "manifest.json")

        if auto_prepare and not os.path.exists(manifest_path):
            logger.info("Staging directory empty -- running preprocessor")
            prep = Preprocessor(config_path)
            prep.run()

        self.engine = DuckDBEngine(
            staging_dir=staging_dir,
            cache_dir=cache_dir,
            session_name=session_name,
        )

        self._score_views: list[str] = []
        self._score_cols: list[str] = []
        self._eval_views: list[str] = []
        self._eval_aliases: list[str] = []
        self._filter_views: list[str] = []
        self._filter_enabled: dict[str, bool] = {}
        self._join_type: str = self.registry.pipeline.join_type
        self._history: list[Action] = []

        self._try_restore_session()

    def _try_restore_session(self) -> None:
        """If the DuckDB file already has a session view, restore state from metadata."""
        try:
            meta_scores = self.engine.get_session_meta("score_views")
            if meta_scores:
                self._score_views = json.loads(meta_scores)
            meta_score_cols = self.engine.get_session_meta("score_cols")
            if meta_score_cols:
                self._score_cols = json.loads(meta_score_cols)
            meta_evals = self.engine.get_session_meta("eval_views")
            if meta_evals:
                self._eval_views = json.loads(meta_evals)
            meta_eval_aliases = self.engine.get_session_meta("eval_aliases")
            if meta_eval_aliases:
                self._eval_aliases = json.loads(meta_eval_aliases)
            meta_filters = self.engine.get_session_meta("filter_views")
            if meta_filters:
                self._filter_views = json.loads(meta_filters)
            meta_join = self.engine.get_session_meta("join_type")
            if meta_join:
                self._join_type = meta_join
        except Exception:
            logger.debug("No prior session state to restore")

    def _persist_session_state(self) -> None:
        """Write current session state to DuckDB metadata tables."""
        self.engine.set_session_meta("score_views", json.dumps(self._score_views))
        self.engine.set_session_meta("score_cols", json.dumps(self._score_cols))
        self.engine.set_session_meta("eval_views", json.dumps(self._eval_views))
        self.engine.set_session_meta("eval_aliases", json.dumps(self._eval_aliases))
        self.engine.set_session_meta("filter_views", json.dumps(self._filter_views))
        self.engine.set_session_meta("join_type", self._join_type)

    def _rebuild_session_view(self) -> None:
        """Rebuild the Layer 2 session view from current state."""
        if not self._score_views:
            return
        active_filters = [
            fv for fv in self._filter_views
            if self._filter_enabled.get(fv, True)
        ]
        self.engine.create_session_view(
            score_views=self._score_views,
            join_type=self._join_type,
            eval_views=self._eval_views if self._eval_views else None,
            filter_views=active_filters if active_filters else None,
        )
        self._persist_session_state()

    # ------------------------------------------------------------------ #
    # Core Operations
    # ------------------------------------------------------------------ #

    def load_base(
        self,
        score_names: list[str],
        join_type: str = "inner",
    ) -> None:
        """Register and join selected prediction tables.

        ``score_names`` can be score column names (looked up via the
        registry) or direct view names already registered.
        """
        start = time.perf_counter()
        self._join_type = join_type
        self._score_views = []
        self._score_cols = []

        for name in score_names:
            entry = self.registry.get_prediction_by_score(name)
            view_name = self.engine.register_prediction(entry)
            if view_name not in self._score_views:
                self._score_views.append(view_name)
            for sc in entry.score_columns:
                if sc not in self._score_cols:
                    self._score_cols.append(sc)

        self._rebuild_session_view()

        elapsed = time.perf_counter() - start
        action = Action(
            type="load_base",
            params={"score_names": score_names, "join_type": join_type},
        )
        self._history.append(action)
        logger.info("Loaded base with %d scores in %.1fs", len(score_names), elapsed)

    def add_score(self, score_name: str, join_type: str | None = None) -> None:
        """Add a new score to the session."""
        entry = self.registry.get_prediction_by_score(score_name)
        view_name = self.engine.register_prediction(entry)

        if view_name in self._score_views:
            logger.warning("Score view %s already in session", view_name)
            return

        self._score_views.append(view_name)
        for sc in entry.score_columns:
            if sc not in self._score_cols:
                self._score_cols.append(sc)

        if join_type:
            self._join_type = join_type

        self._rebuild_session_view()
        self._history.append(Action(
            type="add_score",
            params={"score_name": score_name, "view": view_name},
        ))
        logger.info("Added score: %s", score_name)

    def swap_score(self, old_score: str, new_score: str) -> None:
        """Swap one score for another by redefining the session view."""
        old_entry = self.registry.get_prediction_by_score(old_score)
        new_entry = self.registry.get_prediction_by_score(new_score)

        from ui.engine import _sanitize_name
        old_view = f"score_{_sanitize_name(old_entry.path)}"
        new_view = self.engine.register_prediction(new_entry)

        if old_view in self._score_views:
            idx = self._score_views.index(old_view)
            self._score_views[idx] = new_view
        else:
            self._score_views.append(new_view)

        for sc in old_entry.score_columns:
            if sc in self._score_cols:
                self._score_cols.remove(sc)
        for sc in new_entry.score_columns:
            if sc not in self._score_cols:
                self._score_cols.append(sc)

        self._rebuild_session_view()
        self._history.append(Action(
            type="swap_score",
            params={"old": old_score, "new": new_score},
        ))
        logger.info("Swapped %s -> %s", old_score, new_score)

    def remove_score(self, score_name: str) -> None:
        """Remove a score from the session view."""
        entry = self.registry.get_prediction_by_score(score_name)
        from ui.engine import _sanitize_name
        view_name = f"score_{_sanitize_name(entry.path)}"

        if view_name in self._score_views:
            self._score_views.remove(view_name)
        for sc in entry.score_columns:
            if sc in self._score_cols:
                self._score_cols.remove(sc)

        if self._score_views:
            self._rebuild_session_view()
        self._history.append(Action(
            type="remove_score",
            params={"score_name": score_name, "view": view_name},
        ))
        logger.info("Removed score: %s", score_name)

    def attach_eval(self, alias: str) -> None:
        """Attach an evaluation table to the session via LEFT JOIN."""
        evals = self.registry.list_evaluations()
        matched = None
        for e in evals:
            stem = os.path.basename(e.path).split(".")[0]
            if alias in stem or stem.startswith(alias):
                matched = e
                break
        if matched is None:
            raise KeyError(f"No evaluation table matching alias: {alias}")

        view_name = self.engine.register_evaluation(matched, alias)
        if view_name not in self._eval_views:
            self._eval_views.append(view_name)
        if alias not in self._eval_aliases:
            self._eval_aliases.append(alias)

        self._rebuild_session_view()
        self._history.append(Action(
            type="attach_eval",
            params={"alias": alias, "path": matched.path},
        ))
        logger.info("Attached eval: %s", alias)

    def detach_eval(self, alias: str) -> None:
        """Remove an evaluation table from the session."""
        view_name = f"eval_{alias}"
        if view_name in self._eval_views:
            self._eval_views.remove(view_name)
        if alias in self._eval_aliases:
            self._eval_aliases.remove(alias)

        self._rebuild_session_view()
        self._history.append(Action(
            type="detach_eval",
            params={"alias": alias},
        ))
        logger.info("Detached eval: %s", alias)

    def toggle_filter(self, name: str, enabled: bool) -> None:
        """Enable/disable a filter. Filters are predicates, not row removal.

        The filter column remains in the session view regardless; toggling
        determines whether it's active for eval queries.
        """
        view_name = f"filter_{name}"
        if view_name not in self._filter_views:
            for path in self.registry.list_filters():
                from ui.prepare import _derive_filter_name
                fname = _derive_filter_name(path)
                if fname == name or name in fname:
                    registered = self.engine.register_filter(path, name)
                    if registered:
                        self._filter_views.append(registered)
                    break
            else:
                logger.warning("Filter not found: %s", name)
                return

        self._filter_enabled[view_name] = enabled
        self._rebuild_session_view()
        self._history.append(Action(
            type="toggle_filter",
            params={"name": name, "enabled": enabled},
        ))

    # ------------------------------------------------------------------ #
    # Analysis
    # ------------------------------------------------------------------ #

    def run_stat(
        self,
        stat: str = "enrichment",
        eval_set: str | None = None,
        thresholds: list[float] | None = None,
        **kwargs: Any,
    ) -> tuple[pl.DataFrame, list[dict[str, Any]], pl.DataFrame]:
        """Run eval statistics on the current session data.

        Uses the Arrow bridge to convert session data to a Polars
        LazyFrame, then calls ``eval.run_from_frame()``.
        """
        from biostat_cli.cli import run_from_frame

        if not self.engine.has_cache():
            logger.info("Materializing cache before stat computation...")
            self.engine.materialize_cache()

        lf = self.engine.to_lazy()

        eval_cols: list[str] | None = None
        if eval_set:
            eval_cols = [f"is_pos_{eval_set}"]

        filter_dict: dict[str, str] | None = None
        active_filter_cols = [
            fv for fv, enabled in self._filter_enabled.items() if enabled
        ]
        if active_filter_cols:
            filter_dict = {fv: fv for fv in active_filter_cols}

        results, timings, vsm_cmp = run_from_frame(
            lf,
            eval_level="variant",
            score_cols=self._score_cols,
            stat=stat,
            eval_set=eval_cols,
            filters=filter_dict,
            thresholds=thresholds or self.registry.pipeline.percentile_thresholds,
            **kwargs,
        )
        self._history.append(Action(
            type="run_stat",
            params={"stat": stat, "eval_set": eval_set},
        ))
        return results, timings, vsm_cmp

    def compare(
        self,
        score_cols: list[str] | None = None,
        stat: str = "enrichment",
        eval_set: str | None = None,
        thresholds: list[float] | None = None,
    ) -> pl.DataFrame:
        """Run stat for each score, assemble a comparative DataFrame."""
        cols = score_cols or self._score_cols
        from biostat_cli.cli import run_from_frame

        if not self.engine.has_cache():
            self.engine.materialize_cache()

        lf = self.engine.to_lazy()

        eval_cols_list: list[str] | None = None
        if eval_set:
            eval_cols_list = [f"is_pos_{eval_set}"]

        results, _, _ = run_from_frame(
            lf,
            eval_level="variant",
            score_cols=cols,
            stat=stat,
            eval_set=eval_cols_list,
            thresholds=thresholds or self.registry.pipeline.percentile_thresholds,
        )
        return results

    # ------------------------------------------------------------------ #
    # Persistence & Export
    # ------------------------------------------------------------------ #

    def save(self, path: str | None = None, include_percentiles: bool = True) -> None:
        """Export session table to parquet.

        Optionally adds percentile columns via DuckDB ``PERCENT_RANK()``
        before writing.
        """
        output = path or self.registry.output
        pct_cols = self._score_cols if include_percentiles else None
        self.engine.export_parquet(output, percentile_cols=pct_cols)
        self._history.append(Action(
            type="save",
            params={"path": output, "include_percentiles": include_percentiles},
        ))
        logger.info("Saved session to %s", output)

    def rebuild(self) -> None:
        """Full re-join from scratch. Invalidates cache and rebuilds session view."""
        self.engine.invalidate_cache()
        self._rebuild_session_view()
        self.engine.materialize_cache()
        self._history.append(Action(type="rebuild", params={}))
        logger.info("Rebuilt session from scratch")

    def export_cli_command(self) -> str:
        """Reconstruct an equivalent batch CLI command from session history."""
        pred_paths = []
        for entry in self.registry.list_predictions():
            from ui.engine import _sanitize_name
            view_name = f"score_{_sanitize_name(entry.path)}"
            if view_name in self._score_views:
                pred_paths.append(entry.path)

        parts = [
            "python -m merge.create_vsm_table",
            f"  --prediction_tables {','.join(pred_paths)}",
            f"  --join_type {self._join_type}",
        ]

        eval_paths = []
        for entry in self.registry.list_evaluations():
            stem = os.path.basename(entry.path).split(".")[0]
            for alias in self._eval_aliases:
                if alias in stem:
                    eval_paths.append(entry.path)
                    break

        if eval_paths:
            parts.append(f"  --evaluation_tables {','.join(eval_paths)}")

        linker = self.registry.linker_table
        if linker:
            parts.append(f"  --linker {linker}")

        parts.append(f"  --output {self.registry.output}")
        return " \\\n".join(parts)

    def export_resources_json(self) -> dict[str, Any]:
        """Produce a config dict compatible with eval's resources.json format.

        Includes the merge parameters that produced the current table,
        suitable for use with the batch eval pipeline.
        """
        score_cols = list(self._score_cols)
        eval_entries = []
        for alias in self._eval_aliases:
            for e in self.registry.list_evaluations():
                stem = os.path.basename(e.path).split(".")[0]
                if alias in stem:
                    eval_entries.append({"path": e.path, "alias": alias})
                    break

        filter_info = {}
        for fv in self._filter_views:
            filter_info[fv] = fv

        table_info = {
            "session_export": {
                "Path": self.registry.output,
                "Level": "variant",
                "Score_cols": score_cols,
                "Filters": filter_info,
                "evals": [f"is_pos_{a}" for a in self._eval_aliases],
            }
        }
        return {
            "Table_info": table_info,
            "_session_metadata": {
                "session_name": self.session_name,
                "join_type": self._join_type,
                "score_views": self._score_views,
                "eval_aliases": self._eval_aliases,
                "pipeline": {
                    "join_type": self.registry.pipeline.join_type,
                    "percentile_order": self.registry.pipeline.percentile_order,
                    "reference_score": self.registry.pipeline.reference_score,
                },
            },
        }

    def export_session_log(self) -> str:
        """Produce a human-readable markdown narrative of the session.

        Suitable for lab notebooks or paper methods sections.
        """
        lines = [
            f"# Session Log: {self.session_name}",
            "",
            f"**Config**: `{self.config_path}`",
            f"**Join type**: {self._join_type}",
            f"**Active scores**: {', '.join(self._score_cols) or '(none)'}",
            f"**Attached evals**: {', '.join(self._eval_aliases) or '(none)'}",
            "",
            "## Action History",
            "",
        ]
        for action in self._history:
            ts = action.timestamp.strftime("%Y-%m-%d %H:%M:%S UTC")
            params_str = ", ".join(
                f"{k}={v}" for k, v in action.params.items() if v
            )
            lines.append(f"- **{ts}** `{action.type}`: {params_str}")

        lines.append("")
        lines.append("## Equivalent CLI Command")
        lines.append("")
        lines.append("```bash")
        lines.append(self.export_cli_command())
        lines.append("```")
        return "\n".join(lines)

    def undo(self) -> None:
        """Pop last action from history and reverse it.

        Only supports reversing add/remove/swap/attach/detach operations.
        For load_base, this is a no-op (cannot un-load).
        """
        if not self._history:
            logger.warning("Nothing to undo")
            return

        last = self._history.pop()

        if last.type == "add_score":
            score_name = last.params.get("score_name", "")
            if score_name:
                self.remove_score(score_name)
                self._history.pop()
        elif last.type == "remove_score":
            score_name = last.params.get("score_name", "")
            if score_name:
                self.add_score(score_name)
                self._history.pop()
        elif last.type == "swap_score":
            old = last.params.get("old", "")
            new = last.params.get("new", "")
            if old and new:
                self.swap_score(new, old)
                self._history.pop()
        elif last.type == "attach_eval":
            alias = last.params.get("alias", "")
            if alias:
                self.detach_eval(alias)
                self._history.pop()
        elif last.type == "detach_eval":
            alias = last.params.get("alias", "")
            if alias:
                self.attach_eval(alias)
                self._history.pop()
        else:
            logger.warning("Cannot undo action type: %s", last.type)

    # ------------------------------------------------------------------ #
    # Properties
    # ------------------------------------------------------------------ #

    @property
    def score_cols(self) -> list[str]:
        """Currently active score columns."""
        return list(self._score_cols)

    @property
    def eval_cols(self) -> list[str]:
        """Currently attached eval column names (is_pos_*)."""
        return [f"is_pos_{a}" for a in self._eval_aliases]

    @property
    def history(self) -> list[Action]:
        """Session action history (in-memory)."""
        return list(self._history)

    @property
    def info(self) -> dict[str, Any]:
        """Summary of current session state."""
        return {
            "session_name": self.session_name,
            "score_views": self._score_views,
            "score_cols": self._score_cols,
            "eval_aliases": self._eval_aliases,
            "filter_views": self._filter_views,
            "join_type": self._join_type,
            "has_cache": self.engine.has_cache(),
            "table_info": self.engine.table_info(),
        }

    def close(self) -> None:
        """Close the DuckDB connection."""
        self._persist_session_state()
        self.engine.close()
