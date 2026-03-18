from __future__ import annotations

import sys
from functools import reduce

import polars as pl


JOIN_KEYS = ["chrom", "pos", "ref", "alt"]


def merge_tables(
    tables: list[pl.LazyFrame],
    join_type: str = "inner",
) -> pl.LazyFrame:
    """
    Sequentially join *tables* on genomic keys (chrom, pos, ref, alt).

    Uses a pairwise reduce -- safe for both inner and outer joins when
    each table's key set is unique per row.
    """
    if not tables:
        raise ValueError("No tables provided for merging.")
    if len(tables) == 1:
        return tables[0]

    def _join_pair(left: pl.LazyFrame, right: pl.LazyFrame) -> pl.LazyFrame:
        return left.join(
            right,
            on=JOIN_KEYS,
            how=join_type,
            coalesce=True,
        )

    return reduce(_join_pair, tables)

def merge_tables_pairwise(
    tables: list[pl.LazyFrame],
    score_cols: list[str],
    anchor: str,
) -> tuple[pl.LazyFrame, list[str]]:
    """
    Pairwise merge strategy for prediction score tables.

    Strategy: outer-join all tables once, then create pairwise anchor
    columns via null-masking expressions.  For each non-anchor column *C*,
    ``{anchor}_pairwise_{C}`` contains the anchor value only on rows where
    *C* is non-null (semantically equivalent to an independent inner join
    of the anchor with *C*, but without extra join operations).

    Parameters
    ----------
    tables : list[pl.LazyFrame]
        Prediction table LazyFrames containing keys + score columns.
    score_cols : list[str]
        Names of all score columns across the tables.
    anchor : str
        The score column to use as the anchor for pairwise joining.

    Returns
    -------
    merged : pl.LazyFrame
        The outer-joined result containing keys, the original anchor column,
        each non-anchor column, and each ``{anchor}_pairwise_{C}`` column.
    updated_score_cols : list[str]
        Original *score_cols* plus the new pairwise column names.
    """
    if not tables:
        raise ValueError("No tables provided for merging.")

    pruned_tables = []
    for table in tables:
        schema = table.collect_schema()
        cols = [c for c in schema.names() if c in score_cols]
        if cols:
            pruned_tables.append(table.select(JOIN_KEYS + cols))

    if not pruned_tables:
        raise ValueError("No score columns found in the provided tables.")

    merged = merge_tables(pruned_tables, join_type="outer")

    non_anchor_cols = [c for c in score_cols if c != anchor]

    pairwise_exprs = []
    for c in non_anchor_cols:
        both_non_null = pl.col(c).is_not_null() & pl.col(anchor).is_not_null()
        pairwise_exprs.append(
            pl.when(both_non_null).then(pl.col(anchor)).alias(f"{anchor}_pairwise_{c}")
        )
        pairwise_exprs.append(
            pl.when(both_non_null).then(pl.col(c)).alias(c)
        )

    merged = merged.with_columns(pairwise_exprs)

    pairwise_cols = [f"{anchor}_pairwise_{c}" for c in non_anchor_cols]
    updated_score_cols = list(score_cols) + pairwise_cols

    return merged, updated_score_cols


def aggregate_by_gene(
    lf: pl.LazyFrame,
    score_cols: list[str],
) -> tuple[pl.LazyFrame, list[str]]:
    """
    Aggregate variant-level scores to gene level.

    Groups *lf* by ``ENSG`` and computes mean and max for each column
    in *score_cols*.

    Returns
    -------
    aggregated : pl.LazyFrame
        One row per ENSG with ``{col}_mean`` and ``{col}_max`` columns.
    agg_score_cols : list[str]
        The new aggregated column names.
    """
    agg_exprs: list[pl.Expr] = []
    agg_score_cols: list[str] = []

    for col in score_cols:
        agg_exprs.append(pl.col(col).mean().alias(f"{col}_mean"))
        agg_exprs.append(pl.col(col).max().alias(f"{col}_max"))
        agg_score_cols.append(f"{col}_mean")
        agg_score_cols.append(f"{col}_max")

    aggregated = lf.group_by("ensg").agg(agg_exprs)

    print(
        f"  Aggregated to gene level: {len(agg_score_cols)} agg columns",
        file=sys.stderr,
    )

    return aggregated, agg_score_cols
