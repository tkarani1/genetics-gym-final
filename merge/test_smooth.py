"""
Tests for smooth.py.

Unit tests use synthetic data with mocked structure loaders and require no
reference data. Integration tests require the sir-reference-data directory;
they are skipped automatically if SMOOTH_REFERENCE_DIR is not set.
"""
from __future__ import annotations

import os
from unittest.mock import patch

import numpy as np
import polars as pl
import pytest

from smooth import add_smoothed_columns

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

REFERENCE_DIR = os.environ.get(
    "SMOOTH_REFERENCE_DIR",
    os.path.expanduser(
        "~/Partners HealthCare Dropbox/Hilary Finucane/"
        "Hilary Finucane's files/Home/Research/sir/sir-reference-data"
    ),
)
needs_reference = pytest.mark.skipif(
    not os.path.isdir(REFERENCE_DIR),
    reason="SMOOTH_REFERENCE_DIR not available",
)

SCORE_COLS = ["score_a", "score_b"]


def _make_lf(rows: list[dict]) -> pl.LazyFrame:
    return pl.DataFrame(rows).lazy()


def _synthetic_lf(uniprot_id: str, aa_positions: list[int], scores_a: list[float], scores_b: list[float]) -> pl.LazyFrame:
    """Build a minimal LazyFrame with the columns add_smoothed_columns expects."""
    n = len(aa_positions)
    rows = [
        {
            "chrom": "chr1",
            "pos": 1000 + i,
            "ref": "A",
            "alt": "T",
            "_uniprot_id": uniprot_id,
            "_aa_pos": aa_positions[i],
            "score_a": scores_a[i],
            "score_b": scores_b[i],
        }
        for i in range(n)
    ]
    return pl.DataFrame(rows).lazy()


def _mock_ref(lf: pl.LazyFrame, uniprot_id: str):
    """Return a fake _load_ref_for_chrom result matching the LazyFrame."""
    import pandas as pd
    df = lf.collect().to_pandas()
    ref = pd.DataFrame({
        "_variant_key": "chr1-" + df["pos"].astype(str) + "-" + df["ref"] + "-" + df["alt"],
        "uniprot_id": df["_uniprot_id"],
        "aa_pos": df["_aa_pos"],
    })
    return ref


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

class TestKernelMath:
    """Gaussian kernel weights and weighted averages are computed correctly."""

    def _run(self, dist, pae, scores_a, scores_b, sigma=10.0, pae_cutoff=15.0):
        """Run smoothing with fully mocked structure loaders."""
        n = len(scores_a)
        aa_positions = list(range(1, n + 1))
        uniprot_id = "TESTPROT"
        lf = _synthetic_lf(uniprot_id, aa_positions, scores_a, scores_b)

        def fake_ref(h5_path, chrom, pos_filter):
            return _mock_ref(lf, uniprot_id)

        with (
            patch("smooth._load_ref_for_chrom", side_effect=fake_ref),
            patch("smooth._get_distance_matrix_structure", return_value=dist),
            patch("smooth._get_pae_matrix_structure", return_value=pae),
        ):
            result = add_smoothed_columns(
                lf.select(["chrom", "pos", "ref", "alt", "score_a", "score_b"]),
                SCORE_COLS,
                reference_dir="/fake",
                sigma=sigma,
                pae_cutoff=pae_cutoff,
            ).collect().to_pandas()

        return result

    def test_self_smoothing_tiny_sigma(self):
        """With sigma → 0, each variant should only see itself."""
        dist = np.array([[0.0, 5.0], [5.0, 0.0]])
        pae = np.zeros((2, 2))
        scores_a = [0.1, 0.9]
        scores_b = [0.3, 0.7]

        result = self._run(dist, pae, scores_a, scores_b, sigma=0.001)

        assert result["score_a_smoothed"].iloc[0] == pytest.approx(0.1, abs=1e-4)
        assert result["score_a_smoothed"].iloc[1] == pytest.approx(0.9, abs=1e-4)
        assert result["score_b_smoothed"].iloc[0] == pytest.approx(0.3, abs=1e-4)
        assert result["score_b_smoothed"].iloc[1] == pytest.approx(0.7, abs=1e-4)

    def test_large_sigma_averages(self):
        """With sigma → ∞, all variants converge to the simple mean."""
        dist = np.array([[0.0, 5.0, 10.0],
                         [5.0, 0.0,  5.0],
                         [10.0, 5.0, 0.0]])
        pae = np.zeros((3, 3))
        scores_a = [0.2, 0.4, 0.6]
        scores_b = [0.1, 0.5, 0.9]

        result = self._run(dist, pae, scores_a, scores_b, sigma=1e6)

        expected_a = np.mean(scores_a)
        expected_b = np.mean(scores_b)
        for v in result["score_a_smoothed"]:
            assert v == pytest.approx(expected_a, abs=1e-4)
        for v in result["score_b_smoothed"]:
            assert v == pytest.approx(expected_b, abs=1e-4)

    def test_weighted_average_by_hand(self):
        """Manually verify Gaussian-weighted average for a 3-variant case."""
        # v1 at pos 1, v2 at pos 2 (close), v3 at pos 3 (far from v1)
        dist = np.array([[0.0,  2.0, 50.0],
                         [2.0,  0.0, 50.0],
                         [50.0, 50.0, 0.0]])
        pae = np.zeros((3, 3))
        sigma = 10.0
        scores_a = [0.1, 0.9, 0.5]

        result = self._run(dist, pae, scores_a, [0.0, 0.0, 0.0], sigma=sigma)

        # Compute expected for v1 (aa_pos=1)
        W = np.exp(-((dist / sigma) ** 2))
        idx = np.array([0, 1, 2])
        W_sub = W[np.ix_(idx, idx)]
        expected = (W_sub @ np.array(scores_a)) / W_sub.sum(axis=1)

        for i in range(3):
            assert result["score_a_smoothed"].iloc[i] == pytest.approx(expected[i], abs=1e-6)


class TestPAEMasking:
    """PAE > pae_cutoff sets the corresponding distance to infinity."""

    def test_high_pae_gives_zero_weight(self):
        """A neighbour with PAE > cutoff should not influence the smoothed value."""
        # Two variants: v1 and v2 are close (d=2 Å) but PAE is high → v2 masked
        dist = np.array([[0.0, 2.0], [2.0, 0.0]])
        pae  = np.array([[0.0, 20.0], [20.0, 0.0]])  # PAE > 15 between them
        scores_a = [0.1, 0.9]
        uniprot_id = "TESTPROT"
        n = 2
        lf = _synthetic_lf(uniprot_id, [1, 2], scores_a, [0.0, 0.0])

        def fake_ref(h5_path, chrom, pos_filter):
            return _mock_ref(lf, uniprot_id)

        with (
            patch("smooth._load_ref_for_chrom", side_effect=fake_ref),
            patch("smooth._get_distance_matrix_structure", return_value=dist),
            patch("smooth._get_pae_matrix_structure", return_value=pae),
        ):
            result = add_smoothed_columns(
                lf.select(["chrom", "pos", "ref", "alt", "score_a", "score_b"]),
                SCORE_COLS,
                reference_dir="/fake",
                sigma=10.0,
                pae_cutoff=15.0,
            ).collect().to_pandas()

        # Each variant sees only itself (PAE masks the neighbour)
        assert result["score_a_smoothed"].iloc[0] == pytest.approx(0.1, abs=1e-6)
        assert result["score_a_smoothed"].iloc[1] == pytest.approx(0.9, abs=1e-6)

    def test_low_pae_passes(self):
        """A neighbour with PAE <= cutoff contributes normally."""
        dist = np.array([[0.0, 2.0], [2.0, 0.0]])
        pae  = np.array([[0.0, 5.0], [5.0, 0.0]])  # PAE < 15, no masking
        scores_a = [0.1, 0.9]
        uniprot_id = "TESTPROT"
        lf = _synthetic_lf(uniprot_id, [1, 2], scores_a, [0.0, 0.0])

        def fake_ref(h5_path, chrom, pos_filter):
            return _mock_ref(lf, uniprot_id)

        with (
            patch("smooth._load_ref_for_chrom", side_effect=fake_ref),
            patch("smooth._get_distance_matrix_structure", return_value=dist),
            patch("smooth._get_pae_matrix_structure", return_value=pae),
        ):
            result = add_smoothed_columns(
                lf.select(["chrom", "pos", "ref", "alt", "score_a", "score_b"]),
                SCORE_COLS,
                reference_dir="/fake",
                sigma=10.0,
            ).collect().to_pandas()

        w = np.exp(-((2.0 / 10.0) ** 2))
        expected_v1 = (1.0 * 0.1 + w * 0.9) / (1.0 + w)
        assert result["score_a_smoothed"].iloc[0] == pytest.approx(expected_v1, abs=1e-6)


class TestMissingData:
    """Proteins with any null score yield all-null smoothed columns."""

    def test_null_score_produces_null_smoothed(self):
        uniprot_id = "TESTPROT"
        scores_a = [0.5, None]  # one null
        lf = _synthetic_lf(uniprot_id, [1, 2], scores_a, [0.3, 0.7])
        dist = np.array([[0.0, 5.0], [5.0, 0.0]])

        def fake_ref(h5_path, chrom, pos_filter):
            return _mock_ref(lf, uniprot_id)

        with (
            patch("smooth._load_ref_for_chrom", side_effect=fake_ref),
            patch("smooth._get_distance_matrix_structure", return_value=dist),
            patch("smooth._get_pae_matrix_structure", return_value=None),
        ):
            result = add_smoothed_columns(
                lf.select(["chrom", "pos", "ref", "alt", "score_a", "score_b"]),
                SCORE_COLS,
                reference_dir="/fake",
            ).collect().to_pandas()

        assert result["score_a_smoothed"].isna().all()
        assert result["score_b_smoothed"].isna().all()

    def test_no_null_scores_all_smoothed(self):
        uniprot_id = "TESTPROT"
        scores_a = [0.3, 0.7]
        lf = _synthetic_lf(uniprot_id, [1, 2], scores_a, [0.2, 0.8])
        dist = np.array([[0.0, 5.0], [5.0, 0.0]])

        def fake_ref(h5_path, chrom, pos_filter):
            return _mock_ref(lf, uniprot_id)

        with (
            patch("smooth._load_ref_for_chrom", side_effect=fake_ref),
            patch("smooth._get_distance_matrix_structure", return_value=dist),
            patch("smooth._get_pae_matrix_structure", return_value=None),
        ):
            result = add_smoothed_columns(
                lf.select(["chrom", "pos", "ref", "alt", "score_a", "score_b"]),
                SCORE_COLS,
                reference_dir="/fake",
            ).collect().to_pandas()

        assert result["score_a_smoothed"].notna().all()
        assert result["score_b_smoothed"].notna().all()


class TestOutputShape:
    """Output has original columns plus smoothed columns, nothing dropped."""

    def test_original_columns_preserved(self):
        uniprot_id = "TESTPROT"
        lf = _synthetic_lf(uniprot_id, [1, 2], [0.3, 0.7], [0.2, 0.8])
        dist = np.array([[0.0, 5.0], [5.0, 0.0]])

        def fake_ref(h5_path, chrom, pos_filter):
            return _mock_ref(lf, uniprot_id)

        input_lf = lf.select(["chrom", "pos", "ref", "alt", "score_a", "score_b"])
        with (
            patch("smooth._load_ref_for_chrom", side_effect=fake_ref),
            patch("smooth._get_distance_matrix_structure", return_value=dist),
            patch("smooth._get_pae_matrix_structure", return_value=None),
        ):
            result = add_smoothed_columns(
                input_lf, SCORE_COLS, reference_dir="/fake",
            ).collect()

        schema_names = result.schema.names()
        assert "chrom" in schema_names
        assert "pos" in schema_names
        assert "ref" in schema_names
        assert "alt" in schema_names
        assert "score_a" in schema_names
        assert "score_b" in schema_names
        assert "score_a_smoothed" in schema_names
        assert "score_b_smoothed" in schema_names
        assert len(result) == 2

    def test_unmapped_variants_get_null_smoothed(self):
        """Variants not in the reference get null smoothed columns."""
        lf = _make_lf([
            {"chrom": "chr1", "pos": 9999999, "ref": "A", "alt": "T",
             "score_a": 0.5, "score_b": 0.5},
        ])

        def fake_ref(h5_path, chrom, pos_filter):
            return None  # no mapping found

        with patch("smooth._load_ref_for_chrom", side_effect=fake_ref):
            result = add_smoothed_columns(
                lf, SCORE_COLS, reference_dir="/fake",
            ).collect().to_pandas()

        assert result["score_a_smoothed"].isna().all()
        assert result["score_b_smoothed"].isna().all()


# ---------------------------------------------------------------------------
# Integration tests (require reference data)
# ---------------------------------------------------------------------------

@needs_reference
class TestIntegration:
    """
    Integration tests using real reference data for protein A0A1B0GV90
    (chr1, 354 variants, 54 unique residues) -- the smallest protein in chr1.
    """

    UNIPROT_ID = "A0A1B0GV90"

    def _build_input_lf(self, n_variants: int | None = None, seed: int = 42):
        """Build a synthetic score table with real genomic coordinates."""
        import h5py
        import hdf5plugin  # noqa: F401
        import pandas as pd

        h5_path = os.path.join(REFERENCE_DIR, "all_missense_variants_gr38.h5")
        with h5py.File(h5_path, "r") as f:
            uid = np.array([x.decode("ascii") for x in f["chr1_uniprot_id"][:].flatten()])
            pos = f["chr1_pos"][:]
            ref_alt = f["chr1_ref_alt"][:]
            mask = uid == self.UNIPROT_ID
            positions = pos[mask, 0]
            refs = np.array([x.decode("ascii") for x in ref_alt[mask, 0].flatten()])
            alts = np.array([x.decode("ascii") for x in ref_alt[mask, 1].flatten()])

        if n_variants is not None:
            rng = np.random.default_rng(seed)
            idx = rng.choice(len(positions), size=min(n_variants, len(positions)), replace=False)
            positions, refs, alts = positions[idx], refs[idx], alts[idx]

        rng = np.random.default_rng(seed)
        df = pl.DataFrame({
            "chrom": ["chr1"] * len(positions),
            "pos":   positions.tolist(),
            "ref":   refs.tolist(),
            "alt":   alts.tolist(),
            "score_a": rng.random(len(positions)).tolist(),
            "score_b": rng.random(len(positions)).tolist(),
        })
        return df.lazy()

    def test_output_shape(self):
        lf = self._build_input_lf()
        result = add_smoothed_columns(lf, SCORE_COLS, REFERENCE_DIR).collect()
        assert "score_a_smoothed" in result.schema.names()
        assert "score_b_smoothed" in result.schema.names()
        assert len(result) == 354

    def test_smoothed_values_in_range(self):
        """Smoothed percentile scores should stay in [0, 1]."""
        lf = self._build_input_lf()
        result = add_smoothed_columns(lf, SCORE_COLS, REFERENCE_DIR).collect().to_pandas()
        smoothed_a = result["score_a_smoothed"].dropna()
        assert len(smoothed_a) > 0, "All smoothed values are null"
        assert smoothed_a.between(0.0, 1.0).all(), "Smoothed values outside [0, 1]"

    def test_smoothing_changes_values(self):
        """Smoothed scores should differ from original (scores are not constant)."""
        lf = self._build_input_lf()
        result = add_smoothed_columns(lf, SCORE_COLS, REFERENCE_DIR).collect().to_pandas()
        orig = result["score_a"].dropna()
        smoothed = result["score_a_smoothed"].dropna()
        assert not np.allclose(orig.values, smoothed.values), \
            "Smoothed scores are identical to originals"

    def test_original_columns_unchanged(self):
        """Original score columns must be identical before and after smoothing."""
        lf = self._build_input_lf()
        original = lf.collect().to_pandas()
        result = add_smoothed_columns(lf, SCORE_COLS, REFERENCE_DIR).collect().to_pandas()
        np.testing.assert_array_equal(original["score_a"].values, result["score_a"].values)
        np.testing.assert_array_equal(original["score_b"].values, result["score_b"].values)

    def test_large_sigma_reduces_variance(self):
        """Larger sigma should reduce variance of smoothed scores."""
        lf = self._build_input_lf()
        result_small = add_smoothed_columns(
            lf, SCORE_COLS, REFERENCE_DIR, sigma=1.0
        ).collect().to_pandas()
        result_large = add_smoothed_columns(
            lf, SCORE_COLS, REFERENCE_DIR, sigma=100.0
        ).collect().to_pandas()

        var_small = result_small["score_a_smoothed"].dropna().var()
        var_large = result_large["score_a_smoothed"].dropna().var()
        assert var_large < var_small, \
            f"Expected larger sigma to reduce variance: {var_large:.4f} vs {var_small:.4f}"
