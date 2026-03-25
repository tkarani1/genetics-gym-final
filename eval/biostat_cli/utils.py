"""Shared utility functions for biostat_cli."""

from __future__ import annotations

import polars as pl


def normalize_chromosome_sort_expr(chrom_col: str) -> pl.Expr:
    """
    Create an expression that normalizes chromosome values to sortable integers.

    Handles both "chr" prefixed and unprefixed chromosome names.
    Maps: 1-22 -> 1-22, X -> 23, Y -> 24, others -> 10000.

    Args:
        chrom_col: Name of the chromosome column.

    Returns:
        Polars expression for chromosome sort order.
    """
    chrom_norm = pl.col(chrom_col).cast(pl.Utf8).str.replace(r"(?i)^chr", "").str.to_uppercase()
    chr_num = chrom_norm.cast(pl.Int64, strict=False)
    return (
        pl.when(chr_num.is_not_null())
        .then(chr_num)
        .when(chrom_norm == "X")
        .then(pl.lit(23))
        .when(chrom_norm == "Y")
        .then(pl.lit(24))
        .otherwise(pl.lit(10_000))
    )


def sort_by_genomic_position(
    df: pl.DataFrame,
    chrom_col: str,
    pos_col: str,
    group_sort_cols: list[str] | None = None,
    additional_sort_exprs: list[pl.Expr] | None = None,
) -> pl.DataFrame:
    """
    Sort a DataFrame by genomic position (chromosome, then position).

    Args:
        df: DataFrame to sort.
        chrom_col: Name of the chromosome column.
        pos_col: Name of the position column.
        group_sort_cols: Optional columns to sort by first.
        additional_sort_exprs: Optional additional sort expressions.

    Returns:
        Sorted DataFrame.
    """
    cols = set(df.columns)
    if chrom_col not in cols or pos_col not in cols:
        return df

    chr_sort_expr = normalize_chromosome_sort_expr(chrom_col).alias("__chr_sort")
    pos_sort_expr = pl.col(pos_col).cast(pl.Int64, strict=False).fill_null(9_999_999_999).alias("__pos_sort")

    sort_cols = list(group_sort_cols or [])
    if additional_sort_exprs:
        for i, expr in enumerate(additional_sort_exprs):
            df = df.with_columns(expr.alias(f"__add_sort_{i}"))
            sort_cols.append(f"__add_sort_{i}")

    df = df.with_columns([chr_sort_expr, pos_sort_expr])
    sort_cols.extend(["__chr_sort", "__pos_sort"])
    df = df.sort(sort_cols)

    drop_cols = ["__chr_sort", "__pos_sort"]
    drop_cols.extend([f"__add_sort_{i}" for i in range(len(additional_sort_exprs or []))])
    return df.drop([c for c in drop_cols if c in df.columns])


def missing_category_sort_expr() -> pl.Expr:
    """
    Create sort expression for missing_category column.

    Orders: all_methods (0) < partial_methods (1) < other (2).
    """
    return (
        pl.when(pl.col("missing_category") == "all_methods")
        .then(pl.lit(0))
        .when(pl.col("missing_category") == "partial_methods")
        .then(pl.lit(1))
        .otherwise(pl.lit(2))
    )


__all__ = [
    "normalize_chromosome_sort_expr",
    "sort_by_genomic_position",
    "missing_category_sort_expr",
]
