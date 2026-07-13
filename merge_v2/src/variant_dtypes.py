#!/usr/bin/env python3
"""Shared variant-key / value dtype policy for the ``merge_v2`` build scripts.

The pipeline can store the variant key + values in one of two type *profiles*,
selected by a single boolean ``compact`` flag that every relevant script exposes
as ``--compact-dtypes`` (default **off**) and the orchestrator forwards run-wide:

================  =========================  ==================================
field             default (compact=False)    compact=True
================  =========================  ==================================
``chrom``         ``VARCHAR`` (``'chr1'``..)  ``UTINYINT`` -- ``chr`` stripped;
                                              1-22 keep their number, ``X``->23,
                                              ``Y``->24, ``M``/``MT``->25. Any
                                              contig that does not map (alt /
                                              unplaced contigs, decoys, ...) ->
                                              ``NULL``, so it drops out on the
                                              non-null key filter / inner joins.
``pos``           ``BIGINT``                  ``UINTEGER`` (uint32; human
                                              positions max ~2.5e8, well under
                                              the ~4.29e9 ceiling).
``ref`` / ``alt`` ``VARCHAR``                 ``UTINYINT`` -- the single-base
                                              ASCII byte as a uint8 (the ``|S1``
                                              byte; ``A``=65, ``C``=67, ``G``=71,
                                              ``T``=84), reversible via ``chr()``
                                              (build an ``S1`` view for display).
                                              **SNV-only**: rows with a
                                              multi-character allele (indel / MNV)
                                              are *dropped* (see
                                              :func:`snv_only_predicate`), so every
                                              value fits in one byte.
``score``         ``FLOAT``                   ``DOUBLE``
``is_pos``        ``BOOLEAN``                 ``BOOLEAN`` (1-byte bool)
================  =========================  ==================================

``is_pos`` is ``BOOLEAN`` in both profiles (only its default was ever a bool);
``chrom``, ``pos``, ``ref``/``alt`` and ``score`` are the fields that re-type.
The SNV drop applies only in the compact profile, where the single ``UTINYINT``
byte cannot represent a multi-character allele.

Because downstream tables inherit the key columns unchanged (``SELECT *``) and
then join on ``(chrom, pos, ref, alt)``, the flag must be applied *consistently*
to every table in a run -- mixing a compact table with a default one would make
those joins silently produce no matches. The orchestrator therefore forwards one
run-wide flag to all key-producing steps (the score / eval builders) **and** to
the steps that re-read raw external keys (the linker + filter joins).

Everything here only *emits DuckDB SQL fragments* (or, for the guard, runs a
single validation query); no table is materialized.
"""

from __future__ import annotations

# Chromosome codes for the non-numeric contigs (1-22 keep their own number).
CHR_X, CHR_Y, CHR_M = 23, 24, 25

# argparse help text, shared so every script documents the flag identically.
FLAG_HELP = (
    "Store the variant key + values in the compact numeric profile "
    "(chrom UTINYINT, pos UINTEGER, ref/alt UTINYINT single-byte ASCII "
    "[non-SNV rows dropped], is_pos BOOLEAN, score DOUBLE). Default: the "
    "historical types (chrom VARCHAR, pos BIGINT, ref/alt VARCHAR, "
    "is_pos BOOLEAN, score FLOAT), keeping all rows."
)


def _strip_chr(varchar_expr: str) -> str:
    """Strip a leading ``chr``/``CHR`` prefix from a VARCHAR chrom expression."""
    return f"regexp_replace({varchar_expr}, '^chr', '', 'i')"


def encode_chrom(varchar_expr: str) -> str:
    """Map a VARCHAR chromosome (``'chr1'``/``'1'``/``'chrX'``/``'chrMT'``) to UTINYINT.

    1-22 keep their number, ``X``->23, ``Y``->24, ``M``/``MT``->25. Any contig
    that does not match becomes ``NULL`` (so it falls out on the non-null key
    filter / inner joins downstream).
    """
    s = _strip_chr(varchar_expr)
    return (
        f"CAST(CASE upper({s}) "
        f"WHEN 'X' THEN {CHR_X} "
        f"WHEN 'Y' THEN {CHR_Y} "
        f"WHEN 'M' THEN {CHR_M} "
        f"WHEN 'MT' THEN {CHR_M} "
        f"ELSE TRY_CAST({s} AS SMALLINT) END AS UTINYINT)"
    )


def encode_pos(expr: str) -> str:
    """Cast a position expression to ``UINTEGER`` (uint32)."""
    return f"TRY_CAST({expr} AS UINTEGER)"


def encode_allele(varchar_expr: str) -> str:
    """Encode a single-base VARCHAR allele as its ASCII byte (``UTINYINT`` / uint8).

    This is the ``|S1`` byte reinterpreted as a uint8 (``A``=65, ``C``=67,
    ``G``=71, ``T``=84, ...), so it round-trips to the base via ``chr()`` -- e.g.
    a downstream ``S1`` view for display. Assumes SNV-only input; drop
    multi-character alleles first with :func:`snv_only_predicate` (a single byte
    cannot represent an indel / MNV).
    """
    return f"CAST(ascii({varchar_expr}) AS UTINYINT)"


def key_expr(canon: str, raw_expr: str, compact: bool) -> str:
    """Wrap a *raw* canonical key expression per the active type profile.

    ``raw_expr`` is the caller's existing expression yielding the **default**
    type (VARCHAR chrom, BIGINT pos, VARCHAR ref/alt). With ``compact=False`` it
    is returned unchanged (byte-for-byte the historical SQL); with
    ``compact=True`` ``chrom``/``pos`` are re-typed and ``ref``/``alt`` are
    encoded to their single ASCII byte (``UTINYINT``). Because the encoded allele
    is an integer, the SNV drop (:func:`snv_only_predicate`) must be applied to
    the **raw VARCHAR** allele, before this encoding.
    """
    if not compact:
        return raw_expr
    if canon == "chrom":
        return encode_chrom(raw_expr)
    if canon == "pos":
        return encode_pos(raw_expr)
    if canon in ("ref", "alt"):
        return encode_allele(raw_expr)
    return raw_expr


def score_type(compact: bool) -> str:
    """Storage type for a numeric score / score-derived column."""
    return "DOUBLE" if compact else "FLOAT"


def bool_field_spec(bool_cast_expr: str, compact: bool) -> tuple[str, str]:
    """``(cast_expr, dedup_agg)`` for an ``is_pos``-style boolean label.

    ``is_pos`` is stored as ``BOOLEAN`` in **both** profiles (a 1-byte bool once
    materialized, deduped with ``bool_or``); ``compact`` is accepted only so the
    call sites stay uniform. Kept as a hook in case the profiles diverge again.
    """
    del compact  # is_pos is BOOLEAN in both profiles
    return (bool_cast_expr, "bool_or")


def snv_only_predicate(
    ref_expr: str, alt_expr: str, compact: bool
) -> str | None:
    """A SQL predicate keeping only SNV rows (single-char ref/alt), or ``None``.

    In the compact profile ``ref``/``alt`` are stored as a single ASCII byte, so
    any row whose allele is longer than one character (indel / MNV) is **dropped**
    up front -- truncating it to one byte would collapse distinct variants onto
    one key. Rows with a NULL allele fail the ``length() = 1`` test and are
    dropped too (they carry an incomplete key and are unusable anyway). Returns
    ``None`` in the default profile, so callers apply no extra filtering.

    ``ref_expr`` / ``alt_expr`` must be the **raw VARCHAR** allele expressions and
    the predicate must be applied *before* :func:`encode_allele` turns them into a
    ``UTINYINT`` byte (``length()`` on the encoded integer would be meaningless).
    """
    if not compact:
        return None
    return f"length({ref_expr}) = 1 AND length({alt_expr}) = 1"


def and_where(*predicates: str | None) -> str:
    """Join non-empty predicates with ``AND`` into a single WHERE body.

    Skips ``None``/empty entries, so an inactive :func:`snv_only_predicate`
    simply drops out. Returns ``TRUE`` if nothing is left (a valid no-op WHERE).
    """
    parts = [p for p in predicates if p]
    return " AND ".join(parts) if parts else "TRUE"
