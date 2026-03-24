"""
Spatial smoothing of VSM scores using 3D protein structure distances.

For each variant, smoothed score = Gaussian-kernel-weighted average of scores
across all other variants in the same protein, where weights are based on
CA-atom pairwise distances from AlphaFold PDB structures.

Distances where PAE > pae_cutoff are treated as infinity (weight = 0).
Proteins with any missing score values are skipped (smoothed columns = null).
"""
from __future__ import annotations

import gzip
import json
import os
import re

import h5py
import hdf5plugin  # noqa: F401 -- registers HDF5 compression codecs
import numpy as np
import pandas as pd
import polars as pl
import sklearn.metrics
from Bio.PDB import PDBParser


# ---------------------------------------------------------------------------
# Structure loading utilities (adapted from structure-informed-rvas/utils.py)
# ---------------------------------------------------------------------------

def _get_pairwise_distances(pdb_file: str, *args) -> np.ndarray:
    """Return CA-atom pairwise distance matrix (Å) for a PDB file.

    Optional positional args i, j are 1-based residue range endpoints.
    """
    parser = PDBParser(QUIET=True)
    if pdb_file.endswith(".gz"):
        with gzip.open(pdb_file, "rt") as handle:
            structure = parser.get_structure("protein", handle)
    else:
        with open(pdb_file) as handle:
            structure = parser.get_structure("protein", handle)

    ca_atoms = [
        residue["CA"].get_coord()
        for model in structure
        for chain in model
        for residue in chain
        if "CA" in residue
    ]

    if len(args) > 1:
        ca_atoms = np.array(ca_atoms[args[0] - 1 : args[1]])
    elif len(args) > 0:
        ca_atoms = np.array(ca_atoms[args[0] - 1 :])
    else:
        ca_atoms = np.array(ca_atoms)

    return sklearn.metrics.pairwise_distances(ca_atoms)


def _get_distance_matrix_structure(
    guide_path: str, pdb_dir: str, uniprot_id: str
) -> np.ndarray | None:
    """Return full-protein CA distance matrix (Å), handling multi-domain proteins."""
    info = pd.read_csv(guide_path, sep="\t")
    pdb_files = info.loc[info.pdb_filename.str.contains(uniprot_id), "pdb_filename"]

    if len(pdb_files) == 0:
        return None

    if len(pdb_files) == 1:
        return _get_pairwise_distances(os.path.join(pdb_dir, pdb_files.iloc[0]))

    # Multi-domain: assemble full matrix using "most central" rule for overlaps
    info = info.iloc[pdb_files.index].copy().reset_index()
    info["startAA"] = info.apply(
        lambda x: int(re.findall(r"\d+", x.pos_covered)[0]), axis=1
    )
    info["endAA"] = info.apply(
        lambda x: int(re.findall(r"\d+", x.pos_covered)[1]), axis=1
    )
    nAA = info.endAA.max()
    info["startAA_next"] = info.startAA.shift(periods=-1, fill_value=nAA)
    info["j"] = np.floor((info.endAA + info.startAA_next) / 2).astype(int)
    info["i"] = info.j.shift(periods=1, fill_value=0) + 1

    dist = np.full((nAA, nAA), fill_value=np.inf)
    cum = 0
    for k in range(info.shape[0]):
        path = os.path.join(pdb_dir, info.pdb_filename.values[k])
        i, j = info.startAA[k], info.endAA[k]
        dist[i - 1 : j, i - 1 : j] = _get_pairwise_distances(path)
    for k in range(info.shape[0]):
        path = os.path.join(pdb_dir, info.pdb_filename.values[k])
        i, j = info.i[k], info.j[k]
        i_f = i - (info.startAA[k] - 1)
        j_f = j - (info.startAA[k] - 1)
        n = j_f - i_f + 1
        dist[cum : cum + n, cum : cum + n] = _get_pairwise_distances(path, i_f, j_f)
        cum += n

    return dist


def _get_paes(pae_file: str, *args) -> np.ndarray:
    """Return PAE matrix for a single PAE JSON file."""
    opener = gzip.open(pae_file, "rt") if pae_file.endswith(".gz") else open(pae_file)
    with opener as f:
        data = json.load(f)
    pae = np.array(data[0]["predicted_aligned_error"])

    if len(args) > 1:
        pae = pae[args[0] - 1 : args[1], args[0] - 1 : args[1]]
    elif len(args) > 0:
        pae = pae[args[0] - 1 :, args[0] - 1 :]
    return pae


def _get_pae_matrix_structure(
    guide_path: str, pae_dir: str, uniprot_id: str
) -> np.ndarray | None:
    """Return full-protein PAE matrix, handling multi-domain proteins.

    Returns None if no PAE data is available for this protein.
    """
    info = pd.read_csv(guide_path, sep="\t")
    pae_files = info.loc[
        info.pae_filename.notnull() & info.pae_filename.str.contains(uniprot_id),
        "pae_filename",
    ]

    if len(pae_files) == 0:
        return None

    if len(pae_files) == 1:
        path = os.path.join(pae_dir, pae_files.iloc[0])
        return _get_paes(path) if os.path.exists(path) else None

    # Multi-domain
    info = info.iloc[pae_files.index].copy().reset_index()
    info["startAA"] = info.apply(
        lambda x: int(re.findall(r"\d+", x.pos_covered)[0]), axis=1
    )
    info["endAA"] = info.apply(
        lambda x: int(re.findall(r"\d+", x.pos_covered)[1]), axis=1
    )
    nAA = info.endAA.max()
    info["startAA_next"] = info.startAA.shift(periods=-1, fill_value=nAA)
    info["j"] = np.floor((info.endAA + info.startAA_next) / 2).astype(int)
    info["i"] = info.j.shift(periods=1, fill_value=0) + 1

    pae = np.full((nAA, nAA), fill_value=np.inf)
    cum = 0
    for k in range(info.shape[0]):
        path = os.path.join(pae_dir, info.pae_filename.values[k])
        i, j = info.startAA[k], info.endAA[k]
        pae[i - 1 : j, i - 1 : j] = _get_paes(path)
    for k in range(info.shape[0]):
        path = os.path.join(pae_dir, info.pae_filename.values[k])
        i, j = info.i[k], info.j[k]
        i_f = i - (info.startAA[k] - 1)
        j_f = j - (info.startAA[k] - 1)
        n = j_f - i_f + 1
        pae[cum : cum + n, cum : cum + n] = _get_paes(path, i_f, j_f)
        cum += n

    return np.minimum(pae, pae.T)


# ---------------------------------------------------------------------------
# Reference mapping (adapted from structure-informed-rvas/read_data.py)
# ---------------------------------------------------------------------------

def _load_ref_for_chrom(
    h5_path: str, chrom: str, pos_filter
) -> pd.DataFrame | None:
    """Load variant→protein mapping from HDF5 reference for one chromosome."""
    with h5py.File(h5_path, "r") as f:
        if f"{chrom}_ref_alt" not in f:
            return None
        ref_alt = f[f"{chrom}_ref_alt"][:]
        pdb_filename = f[f"{chrom}_filename"][:]
        uniprot_id = f[f"{chrom}_uniprot_id"][:]
        positions = f[f"{chrom}_pos"][:]

    pos_set = set(pos_filter)
    subset = np.where(pd.Series(positions[:, 0].flatten()).isin(pos_set))[0]

    df = pd.DataFrame({
        "ref":         ref_alt[subset, 0].flatten(),
        "alt":         ref_alt[subset, 1].flatten(),
        "uniprot_id":  uniprot_id[subset].flatten(),
        "pos":         positions[subset, 0].flatten(),
        "aa_pos":      positions[subset, 1].flatten(),
    })
    for col in ["ref", "alt", "uniprot_id"]:
        df[col] = df[col].str.decode("ascii")

    df["_variant_key"] = (
        chrom + "-" + df["pos"].astype(str) + "-" + df["ref"] + "-" + df["alt"]
    )
    return df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def add_smoothed_columns(
    lf: pl.LazyFrame,
    score_columns: list[str],
    reference_dir: str,
    sigma: float = 10.0,
    pae_cutoff: float = 15.0,
) -> pl.LazyFrame:
    """
    Append a ``{col}_smoothed`` column for each column in *score_columns*.

    Smoothing uses a Gaussian kernel over CA-atom pairwise distances:
        w(i,j) = exp(-(d(i,j) / sigma)^2)
    where distances with PAE > pae_cutoff are treated as infinity (w=0).

    Proteins with any null score value across their variants are skipped and
    their smoothed columns are left null.
    """
    h5_path = os.path.join(reference_dir, "all_missense_variants_gr38.h5")
    guide_path = os.path.join(reference_dir, "pdb_pae_file_pos_guide.tsv")
    pdb_dir = os.path.join(reference_dir, "pdb_files")
    pae_dir = os.path.join(reference_dir, "pae_files")

    df = lf.collect().to_pandas().reset_index(drop=True)

    # Build variant key (chrom-pos-ref-alt) matching the HDF5 reference format
    df["_variant_key"] = (
        df["chrom"].astype(str) + "-"
        + df["pos"].astype(str) + "-"
        + df["ref"].astype(str) + "-"
        + df["alt"].astype(str)
    )

    # Map variants to proteins via HDF5 reference
    ref_frames = []
    for chrom, group in df.groupby("chrom"):
        ref = _load_ref_for_chrom(h5_path, str(chrom), group["pos"])
        if ref is not None:
            ref_frames.append(ref[["_variant_key", "uniprot_id", "aa_pos"]])

    # Join protein info; deduplicate keys so merge stays 1-to-1
    smoothed_cols = [f"{col}_smoothed" for col in score_columns]
    if not ref_frames:
        for sc in smoothed_cols:
            df[sc] = np.nan
        return pl.from_pandas(df.drop(columns=["_variant_key"])).lazy()

    ref_df = (
        pd.concat(ref_frames, ignore_index=True)
        .drop_duplicates(subset=["_variant_key"])
    )

    df_joined = df.merge(ref_df, on="_variant_key", how="left")
    for sc in smoothed_cols:
        df_joined[sc] = np.nan

    # Load pdb_pae_file_pos_guide once (used by both distance and PAE loaders)
    guide_info = pd.read_csv(guide_path, sep="\t")

    proteins_mapped = df_joined.dropna(subset=["uniprot_id"])
    n_skipped_missing = 0
    n_skipped_no_structure = 0
    n_smoothed = 0

    for uniprot_id, group in proteins_mapped.groupby("uniprot_id"):
        # Skip proteins with any missing score value
        if group[score_columns].isnull().any().any():
            n_skipped_missing += 1
            continue

        dist = _get_distance_matrix_structure(guide_path, pdb_dir, uniprot_id)
        if dist is None:
            n_skipped_no_structure += 1
            continue

        pae = _get_pae_matrix_structure(guide_path, pae_dir, uniprot_id)
        if pae is not None:
            dist = dist.copy()
            dist[pae > pae_cutoff] = np.inf

        # Gaussian kernel; exp(-inf) = 0 so PAE-masked pairs contribute nothing
        W = np.exp(-((dist / sigma) ** 2))

        aa_idx = group["aa_pos"].values.astype(int) - 1  # 0-based

        if aa_idx.max() >= dist.shape[0]:
            n_skipped_no_structure += 1
            continue

        # W_sub[i, j] = kernel weight between variant i's residue and variant j's
        W_sub = W[np.ix_(aa_idx, aa_idx)]
        row_sums = W_sub.sum(axis=1)

        for col, sc in zip(score_columns, smoothed_cols):
            scores = group[col].values.astype(float)
            df_joined.loc[group.index, sc] = (W_sub @ scores) / row_sums

        n_smoothed += 1

    print(
        f"  Spatial smoothing: {n_smoothed} proteins smoothed, "
        f"{n_skipped_missing} skipped (missing scores), "
        f"{n_skipped_no_structure} skipped (no structure).",
        file=__import__("sys").stderr,
    )

    output_cols = [c for c in df.columns if c != "_variant_key"] + smoothed_cols
    return pl.from_pandas(df_joined[output_cols]).lazy()
