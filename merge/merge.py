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
    Fully lazy, memory-safe merge of score tables.

    Strategy:
    1. Extract only (JOIN_KEYS + relevant score columns) from each table
    2. Perform a single multi-way inner join across tables
    3. Rename the anchor column to ``{anchor}_anchor`` and create
       ``{anchor}_anchor_with_{C}`` columns for each non-anchor column *C*

    Returns
    -------
    merged : pl.LazyFrame
        LazyFrame with JOIN_KEYS, ``{anchor}_anchor``, each non-anchor
        column, and ``{anchor}_anchor_with_{C}`` pairwise columns.
    updated_score_cols : list[str]
        Renamed anchor + non-anchor columns + pairwise column names.
    """

    if not tables:
        raise ValueError("No tables provided for merging.")

    # Step 1: restrict each table to only relevant columns
    pruned_tables = []
    for table in tables:
        schema = table.collect_schema()
        cols = [c for c in schema.names() if c in score_cols]

        if cols:
            pruned_tables.append(table.select(JOIN_KEYS + cols))

    if not pruned_tables:
        raise ValueError("No score columns found in the provided tables.")

    # Step 2: perform a single join across all tables
    def _join(left: pl.LazyFrame, right: pl.LazyFrame) -> pl.LazyFrame:
        return left.join(right, on=JOIN_KEYS, how="inner")

    merged = reduce(_join, pruned_tables)

    # Step 3: rename anchor and define pairwise columns lazily
    non_anchor_cols = [c for c in score_cols if c != anchor]
    anchor_renamed = f"{anchor}_anchor"

    merged = merged.rename({anchor: anchor_renamed})

    pairwise_exprs = [
        pl.col(anchor_renamed).alias(f"{anchor_renamed}_with_{c}")
        for c in non_anchor_cols
    ]

    merged = merged.with_columns(pairwise_exprs)

    updated_score_cols = (
        [anchor_renamed]
        + non_anchor_cols
        + [f"{anchor_renamed}_with_{c}" for c in non_anchor_cols]
    )

    return merged, updated_score_cols