import polars as pl


def add_percentile_columns(
    lf: pl.LazyFrame,
    score_columns: list[str],
) -> pl.LazyFrame:
    """
    Append a ``{col}_percentile`` column for each column in *score_columns*.

    Uses ``rank(method="max")`` divided by the non-null count so that the
    highest-valued row in each column receives a percentile of exactly 1.0
    (equivalent to SQL ``CUME_DIST()``).  Null score values yield null
    percentile values while preserving the row.
    """
    exprs = [
        (
            pl.col(col).rank(method="max")
            / pl.col(col).count()
        ).alias(f"{col}_percentile")
        for col in score_columns
    ]
    return lf.with_columns(exprs)
