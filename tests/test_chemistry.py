"""Unit tests for titration.chemistry — individual chemistry functions."""

import numpy as np
import pytest

from titration.chemistry import (
    calculate_ionic_strength,
    calculate_log_activity,
    calculate_pk_constants,
    d_h2co3_alk,
    d_hac_alk,
    m_h2o,
    m_hpo4,
    m_nh3,
    per,
    per_co3,
    per_h2co3,
    per_hco3,
)


# ---------------------------------------------------------------------------
# Ionic strength
# ---------------------------------------------------------------------------

class TestIonicStrength:
    """Tests for calculate_ionic_strength(tds, dil)."""

    def test_thesis_defaults(self):
        """TDS=3300, dil=5 gives 0.000025 * (3300/5 - 20) = 0.016."""
        result = calculate_ionic_strength(3300.0, 5.0)
        np.testing.assert_allclose(result, 0.016, rtol=1e-6)

    def test_low_tds(self):
        """TDS=100, dil=1 gives 0.000025 * (100 - 20) = 0.002."""
        result = calculate_ionic_strength(100.0, 1.0)
        np.testing.assert_allclose(result, 0.002, rtol=1e-6)

    def test_zero_ionic_strength(self):
        """TDS/dil == 20 yields exactly zero ionic strength."""
        result = calculate_ionic_strength(100.0, 5.0)
        np.testing.assert_allclose(result, 0.0, atol=1e-10)


# ---------------------------------------------------------------------------
# Activity coefficient
# ---------------------------------------------------------------------------

class TestLogActivity:
    """Tests for calculate_log_activity(mue, ktemp)."""

    def test_sign_negative_for_positive_mue(self):
        """Activity coefficient correction should be negative for positive
        ionic strength."""
        result = calculate_log_activity(0.016, 294.0)
        assert result < 0.0

    def test_zero_ionic_strength_gives_zero(self):
        """At zero ionic strength the log activity correction is zero."""
        result = calculate_log_activity(0.0, 298.0)
        np.testing.assert_allclose(result, 0.0, atol=1e-10)


# ---------------------------------------------------------------------------
# pK constants
# ---------------------------------------------------------------------------

class TestPkConstants:
    """Tests for calculate_pk_constants(temperature, tds, dil)."""

    def test_textbook_pk1_pk2_at_25c(self):
        """At 25 C in near-pure water, pK1 ~ 6.35 and pK2 ~ 10.33."""
        pks = calculate_pk_constants(25.0, 21.0, 1.0)
        np.testing.assert_allclose(pks["pk1"], 6.35, atol=0.01)
        np.testing.assert_allclose(pks["pk2"], 10.33, atol=0.01)

    def test_activity_correction_shifts_pk_values(self, thesis_pk_constants):
        """Activity-corrected pk11/pk22 should differ from raw pk1/pk2."""
        pks = thesis_pk_constants
        assert pks["pk11"] != pks["pk1"]
        assert pks["pk22"] != pks["pk2"]

    def test_logf1_is_negative(self, thesis_pk_constants):
        """logf1 should be negative for positive ionic strength."""
        assert thesis_pk_constants["logf1"] < 0.0

    def test_logf2_is_four_times_logf1(self, thesis_pk_constants):
        """logf2 = 4 * logf1 by definition."""
        pks = thesis_pk_constants
        np.testing.assert_allclose(pks["logf2"], 4.0 * pks["logf1"], rtol=1e-6)

    def test_tds_below_20_clamped_to_21(self):
        """TDS < 20 is clamped to 21 — result must match TDS=21 exactly."""
        pks_15 = calculate_pk_constants(21.0, 15.0, 1.0)
        pks_21 = calculate_pk_constants(21.0, 21.0, 1.0)
        for key in pks_15:
            np.testing.assert_allclose(
                pks_15[key], pks_21[key], rtol=1e-12,
                err_msg=f"Mismatch on {key} for TDS clamping",
            )

    def test_tds_at_threshold_not_clamped(self):
        """TDS=20 is at the boundary (>= 20) and should NOT be clamped."""
        pks_20 = calculate_pk_constants(21.0, 20.0, 1.0)
        pks_21 = calculate_pk_constants(21.0, 21.0, 1.0)
        # They should differ because TDS=20 uses 20, not 21
        assert pks_20["pk11"] != pks_21["pk11"]

    def test_all_twelve_keys_present(self, thesis_pk_constants):
        """The returned dict must contain all expected keys."""
        expected_keys = {
            "pk1", "pk2", "pka", "pkn", "pkp",
            "pk11", "pk22", "pkaa", "pknn", "pkpp",
            "logf1", "logf2",
        }
        assert set(thesis_pk_constants.keys()) == expected_keys


# ---------------------------------------------------------------------------
# Species fractions
# ---------------------------------------------------------------------------

class TestSpeciesFractions:
    """Tests for per_h2co3, per_hco3, per_co3."""

    @pytest.mark.parametrize("ph", [4, 5, 6, 7, 8, 9, 10, 11, 12])
    def test_fractions_sum_to_one(self, thesis_pk_constants, ph):
        """H2CO3 + HCO3 + CO3 must equal 1.0 at every pH."""
        pk11 = thesis_pk_constants["pk11"]
        pk22 = thesis_pk_constants["pk22"]
        total = (
            per_h2co3(ph, pk11, pk22)
            + per_hco3(ph, pk11, pk22)
            + per_co3(ph, pk11, pk22)
        )
        np.testing.assert_allclose(total, 1.0, rtol=1e-12)

    def test_very_low_ph_dominated_by_h2co3(self, thesis_pk_constants):
        """At pH 2, virtually all carbonate is H2CO3*."""
        pk11 = thesis_pk_constants["pk11"]
        pk22 = thesis_pk_constants["pk22"]
        np.testing.assert_allclose(
            per_h2co3(2.0, pk11, pk22), 1.0, atol=1e-4,
        )

    def test_very_high_ph_dominated_by_co3(self, thesis_pk_constants):
        """At pH 13, virtually all carbonate is CO3^2-."""
        pk11 = thesis_pk_constants["pk11"]
        pk22 = thesis_pk_constants["pk22"]
        np.testing.assert_allclose(
            per_co3(13.0, pk11, pk22), 1.0, atol=0.002,
        )

    def test_at_pk11_h2co3_equals_hco3(self, thesis_pk_constants):
        """At pH = pk11, H2CO3* and HCO3- fractions are approximately equal."""
        pk11 = thesis_pk_constants["pk11"]
        pk22 = thesis_pk_constants["pk22"]
        h2co3_frac = per_h2co3(pk11, pk11, pk22)
        hco3_frac = per_hco3(pk11, pk11, pk22)
        np.testing.assert_allclose(h2co3_frac, hco3_frac, rtol=1e-6)


# ---------------------------------------------------------------------------
# Generic weak-acid fraction: per()
# ---------------------------------------------------------------------------

class TestPer:
    """Tests for the generic deprotonated fraction per(ph, pkk)."""

    def test_half_dissociated_at_pka(self):
        """At pH = pKa the fraction must be exactly 0.5."""
        np.testing.assert_allclose(per(5.0, 5.0), 0.5, rtol=1e-12)

    def test_fully_dissociated_at_high_ph(self):
        """At pH >> pKa the fraction approaches 1.0."""
        np.testing.assert_allclose(per(15.0, 5.0), 1.0, atol=1e-9)

    def test_undissociated_at_low_ph(self):
        """At pH << pKa the fraction approaches 0.0."""
        np.testing.assert_allclose(per(1.0, 10.0), 0.0, atol=1e-8)


# ---------------------------------------------------------------------------
# Alkalinity differences
# ---------------------------------------------------------------------------

class TestAlkalinityDifferences:
    """Tests for d_h2co3_alk and d_hac_alk."""

    def test_d_h2co3_alk_same_ph_gives_zero(self, thesis_pk_constants):
        """Same pH in both points yields zero alkalinity difference."""
        pk11 = thesis_pk_constants["pk11"]
        pk22 = thesis_pk_constants["pk22"]
        result = d_h2co3_alk(7.0, 7.0, pk11, pk22)
        np.testing.assert_allclose(result, 0.0, atol=1e-10)

    def test_d_hac_alk_same_ph_gives_zero(self, thesis_pk_constants):
        """Same pH in both points yields zero acetic-acid alkalinity diff."""
        pkaa = thesis_pk_constants["pkaa"]
        result = d_hac_alk(7.0, 7.0, pkaa)
        np.testing.assert_allclose(result, 0.0, atol=1e-10)

    def test_d_h2co3_alk_antisymmetric(self, thesis_pk_constants):
        """Swapping pH arguments negates the result."""
        pk11 = thesis_pk_constants["pk11"]
        pk22 = thesis_pk_constants["pk22"]
        fwd = d_h2co3_alk(7.0, 5.0, pk11, pk22)
        rev = d_h2co3_alk(5.0, 7.0, pk11, pk22)
        np.testing.assert_allclose(fwd, -rev, rtol=1e-12)

    def test_d_hac_alk_antisymmetric(self, thesis_pk_constants):
        """Swapping pH arguments negates the result."""
        pkaa = thesis_pk_constants["pkaa"]
        fwd = d_hac_alk(7.0, 5.0, pkaa)
        rev = d_hac_alk(5.0, 7.0, pkaa)
        np.testing.assert_allclose(fwd, -rev, rtol=1e-12)


# ---------------------------------------------------------------------------
# Mass-balance corrections
# ---------------------------------------------------------------------------

class TestMassBalanceCorrections:
    """Tests for m_nh3, m_hpo4, and m_h2o."""

    def test_m_nh3_zero_nitrogen(self, thesis_pk_constants):
        """With nt=0, nitrogen mass balance returns zero."""
        pknn = thesis_pk_constants["pknn"]
        result = m_nh3(7.0, 6.0, 0.0, 5.0, 50.0, pknn)
        np.testing.assert_allclose(result, 0.0, atol=1e-10)

    def test_m_hpo4_zero_phosphorus(self, thesis_pk_constants):
        """With pt=0, phosphorus mass balance returns zero."""
        pkpp = thesis_pk_constants["pkpp"]
        result = m_hpo4(7.0, 6.0, 0.0, 5.0, 50.0, pkpp)
        np.testing.assert_allclose(result, 0.0, atol=1e-10)

    def test_m_h2o_same_ph_and_volumes_gives_zero(self, thesis_pk_constants):
        """Same pH and same volume for both points gives exactly zero."""
        logf1 = thesis_pk_constants["logf1"]
        result = m_h2o(2.0, 2.0, 7.0, 7.0, 50.0, logf1)
        np.testing.assert_allclose(result, 0.0, atol=1e-10)

    def test_m_nh3_nonzero_nitrogen(self, thesis_pk_constants):
        """With nonzero nt, nitrogen mass balance is nonzero when pH differs."""
        pknn = thesis_pk_constants["pknn"]
        result = m_nh3(7.0, 6.0, 100.0, 5.0, 50.0, pknn)
        assert result != 0.0

    def test_m_hpo4_nonzero_phosphorus(self, thesis_pk_constants):
        """With nonzero pt, phosphorus mass balance is nonzero when pH differs."""
        pkpp = thesis_pk_constants["pkpp"]
        result = m_hpo4(7.0, 6.0, 100.0, 5.0, 50.0, pkpp)
        assert result != 0.0
