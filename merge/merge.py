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
