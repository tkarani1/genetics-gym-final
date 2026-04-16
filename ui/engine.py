"""DuckDB star schema engine with four-layer architecture.

Layer 1: Canonical views -- one VIEW per staged source file (lazy, zero I/O).
Layer 2: Session view   -- composable JOIN of selected Layer 1 views.
Layer 3: Materialized cache -- CREATE TABLE AS on demand.
Layer 4: Export         -- COPY with PERCENT_RANK() window functions.

Assumes clean staged inputs from prepare.py. No schema probing,
key normalization, or format handling.
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Any

import duckdb
import polars as pl

from ui.registry import EvaluationEntry, PredictionEntry

logger = logging.getLogger(__name__)

VARIANT_KEYS = ["chrom", "pos", "ref", "alt"]
GENE_KEY = "ensg"
_SESSION_VIEW = "session_current"
_SESSION_CACHE = "session_cache"
_GENE_VIEW = "session_gene"


def _sanitize_name(raw: str) -> str:
    """Convert a file path or arbitrary string to a valid SQL identifier."""
    stem = Path(raw).stem
    clean = re.sub(r"[^a-zA-Z0-9_]", "_", stem)
    if clean and clean[0].isdigit():
        clean = f"t_{clean}"
    return clean.lower()


class DuckDBEngine:
    """DuckDB star schema engine for session-based data exploration.

    Each engine instance owns a DuckDB connection backed by a persistent
    ``.duckdb`` file in ``cache_dir``. The disk-spill directory is
    configured to ``{cache_dir}/tmp`` so large joins don't OOM.
    """

    def __init__(
        self,
        staging_dir: str = "ui/local_cache/staging",
        cache_dir: str = "ui/local_cache",
        session_name: str = "default",
    ) -> None:
        self.staging_dir = staging_dir
        self.cache_dir = cache_dir
        self.session_name = session_name

        db_path = os.path.join(cache_dir, f"{session_name}.duckdb")
        os.makedirs(cache_dir, exist_ok=True)

        tmp_dir = os.path.join(cache_dir, "tmp")
        os.makedirs(tmp_dir, exist_ok=True)

        self.conn = duckdb.connect(db_path)
        self.conn.execute(f"SET temp_directory = '{tmp_dir}'")

        self._init_meta_tables()
        self._registered_views: dict[str, str] = {}

    def _init_meta_tables(self) -> None:
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS _meta_sources (
                name       VARCHAR PRIMARY KEY,
                path       VARCHAR NOT NULL,
                role       VARCHAR NOT NULL,
                level      VARCHAR NOT NULL,
                columns    VARCHAR[],
                registered TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS _meta_session (
                key   VARCHAR PRIMARY KEY,
                value VARCHAR
            )
        """)
        self.conn.execute("""
            CREATE SEQUENCE IF NOT EXISTS _meta_history_seq START 1
        """)
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS _meta_history (
                seq       INTEGER DEFAULT nextval('_meta_history_seq'),
                action    VARCHAR NOT NULL,
                params    JSON,
                timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

    def _record_source(
        self,
        name: str,
        path: str,
        role: str,
        level: str,
        columns: list[str],
    ) -> None:
        self.conn.execute(
            """
            INSERT OR REPLACE INTO _meta_sources (name, path, role, level, columns)
            VALUES (?, ?, ?, ?, ?)
            """,
            [name, path, role, level, columns],
        )

    def _record_action(self, action: str, params: dict[str, Any] | None = None) -> None:
        import json as json_mod

        self.conn.execute(
            "INSERT INTO _meta_history (action, params) VALUES (?, ?::JSON)",
            [action, json_mod.dumps(params or {})],
        )

    # ------------------------------------------------------------------ #
    # Layer 1: Canonical Views
    # ------------------------------------------------------------------ #

    def register_prediction(self, entry: PredictionEntry) -> str:
        """Register a clean staged prediction file as a DuckDB VIEW."""
        name = f"score_{_sanitize_name(entry.path)}"
        staged = self._resolve_staged_path(entry.path, "predictions")

        keys = ", ".join(VARIANT_KEYS if entry.level == "variant" else [GENE_KEY])
        cols = ", ".join(entry.score_columns)

        self.conn.execute(f"""
            CREATE OR REPLACE VIEW {name} AS
            SELECT {keys}, {cols}
            FROM read_parquet('{staged}')
        """)
        self._registered_views[name] = entry.level
        self._record_source(name, staged, "prediction", entry.level, entry.score_columns)
        logger.info("Registered prediction view: %s (%s)", name, cols)
        return name

    def register_evaluation(self, entry: EvaluationEntry, alias: str) -> str:
        """Register a clean staged eval file as a DuckDB VIEW.

        The label column is aliased as ``is_pos_{alias}`` so multiple
        eval sets can coexist in the session view.
        """
        name = f"eval_{_sanitize_name(alias)}"
        staged = self._resolve_staged_path(entry.path, "evaluations")

        if entry.level == "variant":
            keys = ", ".join(VARIANT_KEYS)
            self.conn.execute(f"""
                CREATE OR REPLACE VIEW {name} AS
                SELECT {keys}, is_pos AS is_pos_{alias}
                FROM read_parquet('{staged}')
            """)
        else:
            self.conn.execute(f"""
                CREATE OR REPLACE VIEW {name} AS
                SELECT * FROM read_parquet('{staged}')
            """)

        self._registered_views[name] = entry.level
        self._record_source(name, staged, "evaluation", entry.level, [f"is_pos_{alias}"])
        logger.info("Registered evaluation view: %s", name)
        return name

    def register_filter(self, path: str, name: str | None = None) -> str:
        """Register a clean staged filter file as a DuckDB VIEW.

        The view adds a boolean ``filter_{name}`` column set to TRUE
        for all rows present in the filter.
        """
        if name is None:
            name = _sanitize_name(path)
        view_name = f"filter_{name}"
        staged = self._resolve_staged_path(path, "filters")

        schema = self.conn.execute(
            f"SELECT * FROM read_parquet('{staged}') LIMIT 0"
        ).description
        col_names = [desc[0] for desc in schema]

        if all(k in col_names for k in VARIANT_KEYS):
            keys = ", ".join(VARIANT_KEYS)
            level = "variant"
        elif GENE_KEY in col_names:
            keys = GENE_KEY
            level = "gene"
        else:
            logger.warning("Cannot determine filter level for %s, skipping", path)
            return ""

        self.conn.execute(f"""
            CREATE OR REPLACE VIEW {view_name} AS
            SELECT DISTINCT {keys}, TRUE AS {view_name}
            FROM read_parquet('{staged}')
        """)
        self._registered_views[view_name] = level
        self._record_source(view_name, staged, "filter", level, [view_name])
        logger.info("Registered filter view: %s", view_name)
        return view_name

    def register_linker(self, path: str | None = None) -> str:
        """Register the linker as a DuckDB VIEW."""
        if path is None:
            staged = os.path.join(self.staging_dir, "linker.parquet")
        else:
            staged = self._resolve_staged_path(path, ".")

        self.conn.execute(f"""
            CREATE OR REPLACE VIEW linker AS
            SELECT * FROM read_parquet('{staged}')
        """)
        self._registered_views["linker"] = "variant"
        self._record_source("linker", staged, "linker", "variant", VARIANT_KEYS + [GENE_KEY])
        logger.info("Registered linker view")
        return "linker"

    def _resolve_staged_path(self, source_path: str, subdir: str) -> str:
        """Find the staged parquet for a source path."""
        basename = Path(source_path).name
        for ext in (".tsv.bgz", ".tsv.gz", ".tsv"):
            if basename.lower().endswith(ext):
                basename = basename[: len(basename) - len(ext)] + ".parquet"
                break
        if not basename.endswith(".parquet"):
            basename = Path(basename).stem + ".parquet"

        candidates = [
            os.path.join(self.staging_dir, subdir, basename),
            os.path.join(self.staging_dir, basename),
        ]
        for c in candidates:
            if os.path.exists(c):
                return c

        return candidates[0]

    # ------------------------------------------------------------------ #
    # Layer 2: Session View
    # ------------------------------------------------------------------ #

    def create_session_view(
        self,
        score_views: list[str],
        join_type: str = "INNER",
        eval_views: list[str] | None = None,
        filter_views: list[str] | None = None,
    ) -> None:
        """Create/replace the session view by joining selected Layer 1 views.

        Score swap = call this again with a different ``score_views`` list.
        """
        if not score_views:
            raise ValueError("At least one score view is required")

        join_kw = join_type.upper()
        if join_kw not in ("INNER", "LEFT", "FULL OUTER"):
            raise ValueError(f"Unsupported join type: {join_type}")

        first_level = self._registered_views.get(score_views[0], "variant")
        if first_level == "variant":
            using_clause = f"USING ({', '.join(VARIANT_KEYS)})"
        else:
            using_clause = f"USING ({GENE_KEY})"

        parts = [f"SELECT * FROM {score_views[0]}"]
        for sv in score_views[1:]:
            parts.append(f"{join_kw} JOIN {sv} {using_clause}")

        if eval_views:
            for ev in eval_views:
                parts.append(f"LEFT JOIN {ev} {using_clause}")

        if filter_views:
            for fv in filter_views:
                parts.append(f"LEFT JOIN {fv} {using_clause}")

        sql = f"CREATE OR REPLACE VIEW {_SESSION_VIEW} AS\n" + "\n".join(parts)
        self.conn.execute(sql)
        self.invalidate_cache()
        self._record_action("create_session_view", {
            "score_views": score_views,
            "join_type": join_type,
            "eval_views": eval_views or [],
            "filter_views": filter_views or [],
        })
        logger.info("Created session view with %d score + %d eval + %d filter views",
                     len(score_views), len(eval_views or []), len(filter_views or []))

    # ------------------------------------------------------------------ #
    # Layer 3: Materialized Cache
    # ------------------------------------------------------------------ #

    def materialize_cache(self) -> None:
        """Materialize the session view into a physical table for fast queries."""
        self.invalidate_cache()
        self.conn.execute(
            f"CREATE TABLE {_SESSION_CACHE} AS SELECT * FROM {_SESSION_VIEW}"
        )
        self._record_action("materialize_cache")
        logger.info("Materialized session cache")

    def invalidate_cache(self) -> None:
        """Drop the materialized cache if it exists."""
        self.conn.execute(f"DROP TABLE IF EXISTS {_SESSION_CACHE}")

    def has_cache(self) -> bool:
        """Check if the materialized cache exists."""
        result = self.conn.execute(
            "SELECT COUNT(*) FROM information_schema.tables "
            f"WHERE table_name = '{_SESSION_CACHE}'"
        ).fetchone()
        return result is not None and result[0] > 0

    # ------------------------------------------------------------------ #
    # Layer 4: Export
    # ------------------------------------------------------------------ #

    def export_parquet(
        self,
        output_path: str,
        percentile_cols: list[str] | None = None,
    ) -> None:
        """Export the session (from cache or view) to a parquet file.

        If ``percentile_cols`` is provided, ``PERCENT_RANK()`` window
        functions are added for each column.
        """
        source = _SESSION_CACHE if self.has_cache() else _SESSION_VIEW

        if percentile_cols:
            pct_exprs = ", ".join(
                f"PERCENT_RANK() OVER (ORDER BY {col}) AS {col}_percentile"
                for col in percentile_cols
            )
            sql = f"SELECT *, {pct_exprs} FROM {source}"
        else:
            sql = f"SELECT * FROM {source}"

        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        self.conn.execute(
            f"COPY ({sql}) TO '{output_path}' (FORMAT PARQUET, COMPRESSION ZSTD)"
        )
        self._record_action("export_parquet", {
            "output_path": output_path,
            "percentile_cols": percentile_cols or [],
        })
        logger.info("Exported to %s", output_path)

    # ------------------------------------------------------------------ #
    # Gene Aggregation
    # ------------------------------------------------------------------ #

    def create_gene_aggregation_view(
        self,
        score_cols: list[str],
        collapse: bool = True,
    ) -> None:
        """Create a gene-level aggregation view from the session data.

        Requires the linker view to be registered. Joins session_current
        with linker via variant keys, then groups by ensg.
        """
        if collapse:
            agg_exprs = ", ".join(
                f"MAX({col}) AS {col}" for col in score_cols
            )
        else:
            agg_exprs = ", ".join(
                f"AVG({col}) AS {col}_mean, MAX({col}) AS {col}_max, "
                f"MIN({col}) AS {col}_min, COUNT({col}) AS {col}_count"
                for col in score_cols
            )

        source = _SESSION_CACHE if self.has_cache() else _SESSION_VIEW

        self.conn.execute(f"""
            CREATE OR REPLACE VIEW {_GENE_VIEW} AS
            SELECT linker.{GENE_KEY}, {agg_exprs}
            FROM {source}
            JOIN linker USING ({', '.join(VARIANT_KEYS)})
            GROUP BY linker.{GENE_KEY}
        """)
        self._record_action("create_gene_aggregation_view", {
            "score_cols": score_cols,
            "collapse": collapse,
        })
        logger.info("Created gene aggregation view (collapse=%s)", collapse)

    # ------------------------------------------------------------------ #
    # Arrow Bridge: DuckDB → Polars
    # ------------------------------------------------------------------ #

    def to_polars(self, columns: list[str] | None = None) -> pl.DataFrame:
        """Read from cache (or view) and return a Polars DataFrame."""
        source = _SESSION_CACHE if self.has_cache() else _SESSION_VIEW
        if columns:
            col_str = ", ".join(columns)
            return self.conn.execute(f"SELECT {col_str} FROM {source}").pl()
        return self.conn.execute(f"SELECT * FROM {source}").pl()

    def to_lazy(self, columns: list[str] | None = None) -> pl.LazyFrame:
        """Same as ``to_polars()`` but returns a ``.lazy()`` for eval consumption."""
        return self.to_polars(columns).lazy()

    def to_polars_gene(self) -> pl.DataFrame:
        """Read from the gene aggregation view."""
        return self.conn.execute(f"SELECT * FROM {_GENE_VIEW}").pl()

    # ------------------------------------------------------------------ #
    # Query Helpers
    # ------------------------------------------------------------------ #

    def execute_sql(self, sql: str) -> pl.DataFrame:
        """Run arbitrary SQL and return as a Polars DataFrame."""
        return self.conn.execute(sql).pl()

    def row_count(self, table: str | None = None) -> int:
        """Row count of session cache/view or a named table."""
        target = table or (_SESSION_CACHE if self.has_cache() else _SESSION_VIEW)
        result = self.conn.execute(f"SELECT COUNT(*) FROM {target}").fetchone()
        return result[0] if result else 0

    def table_info(self) -> dict[str, Any]:
        """Column names, types, and row count of the current session."""
        source = _SESSION_CACHE if self.has_cache() else _SESSION_VIEW
        try:
            desc = self.conn.execute(f"SELECT * FROM {source} LIMIT 0").description
            columns = {d[0]: d[1] for d in desc}
            count = self.row_count()
            return {"columns": columns, "row_count": count, "source": source}
        except duckdb.CatalogException:
            return {"columns": {}, "row_count": 0, "source": None}

    def list_views(self) -> list[str]:
        """List all user-defined views (Layer 1 + session views)."""
        result = self.conn.execute(
            "SELECT table_name FROM information_schema.tables WHERE table_type = 'VIEW'"
        ).fetchall()
        return [r[0] for r in result]

    # ------------------------------------------------------------------ #
    # Session Management
    # ------------------------------------------------------------------ #

    def list_cached_sessions(self) -> list[str]:
        """Enumerate .duckdb files in cache_dir."""
        sessions = []
        for f in os.listdir(self.cache_dir):
            if f.endswith(".duckdb"):
                sessions.append(f.replace(".duckdb", ""))
        return sorted(sessions)

    def get_session_meta(self, key: str) -> str | None:
        """Read a value from the session metadata table."""
        result = self.conn.execute(
            "SELECT value FROM _meta_session WHERE key = ?", [key]
        ).fetchone()
        return result[0] if result else None

    def set_session_meta(self, key: str, value: str) -> None:
        """Write a value to the session metadata table."""
        self.conn.execute(
            "INSERT OR REPLACE INTO _meta_session (key, value) VALUES (?, ?)",
            [key, value],
        )

    def get_history(self, limit: int = 50) -> list[dict[str, Any]]:
        """Return recent session history entries."""
        rows = self.conn.execute(
            f"SELECT seq, action, params, timestamp FROM _meta_history "
            f"ORDER BY seq DESC LIMIT {limit}"
        ).fetchall()
        return [
            {"seq": r[0], "action": r[1], "params": r[2], "timestamp": str(r[3])}
            for r in rows
        ]

    def close(self) -> None:
        """Close the DuckDB connection."""
        self.conn.close()
