"""Functional tests for titration.chemistry — pK constants, species fractions, activity coefficients."""

import numpy as np
import pytest


class TestIonicStrength:
    """Verify calculate_ionic_strength function."""

    def test_thesis_defaults(self):
        from titration.chemistry import calculate_ionic_strength

        # TDS=3300, dil=5.0 (50/10): 0.000025 * (3300/5 - 20) = 0.000025 * 640 = 0.016
        mu = calculate_ionic_strength(3300.0, 5.0)
        np.testing.assert_allclose(mu, 0.016, rtol=1e-10)

    def test_high_tds(self):
        from titration.chemistry import calculate_ionic_strength

        # TDS=10000, dil=5.0: 0.000025 * (10000/5 - 20) = 0.000025 * 1980 = 0.0495
        mu = calculate_ionic_strength(10000.0, 5.0)
        np.testing.assert_allclose(mu, 0.0495, rtol=1e-10)

    def test_undiluted(self):
        from titration.chemistry import calculate_ionic_strength

        # TDS=3300, dil=1.0: 0.000025 * (3300 - 20) = 0.000025 * 3280 = 0.082
        mu = calculate_ionic_strength(3300.0, 1.0)
        np.testing.assert_allclose(mu, 0.082, rtol=1e-10)

    def test_minimum_tds_threshold(self):
        from titration.chemistry import calculate_ionic_strength

        # When TDS is close to dil*20, ionic strength approaches 0
        mu = calculate_ionic_strength(100.0, 5.0)
        # 0.000025 * (100/5 - 20) = 0.000025 * 0 = 0
        np.testing.assert_allclose(mu, 0.0, atol=1e-15)


class TestLogActivity:
    """Verify calculate_log_activity (Debye-Huckel activity coefficient)."""

    def test_thesis_defaults(self):
        from titration.chemistry import calculate_log_activity

        # mu=0.016, ktemp=294.0 (21+273)
        logf = calculate_log_activity(0.016, 294.0)
        # Should be negative for positive ionic strength
        assert logf < 0.0
        np.testing.assert_allclose(logf, -0.05616400115070987, rtol=1e-6)

    def test_zero_ionic_strength(self):
        from titration.chemistry import calculate_log_activity

        # mu=0: sqrt(0)/(1+sqrt(0)) - 0.3*0 = 0
        logf = calculate_log_activity(0.0, 294.0)
        np.testing.assert_allclose(logf, 0.0, atol=1e-15)


class TestPkConstants:
    """Verify calculate_pk_constants returns correct pK values at thesis conditions."""

    def test_thesis_defaults_structure(self):
        from titration.chemistry import calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        expected_keys = {"pk1", "pk2", "pka", "pkn", "pkp",
                         "pk11", "pk22", "pkaa", "pknn", "pkpp",
                         "logf1", "logf2"}
        assert set(pks.keys()) == expected_keys

    def test_pk_values_reasonable_range(self):
        from titration.chemistry import calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        # Carbonate pK1 should be around 6.3-6.4
        assert 6.0 < pks["pk1"] < 6.8
        # Carbonate pK2 should be around 10.2-10.4
        assert 9.8 < pks["pk2"] < 10.8
        # Acetic acid pKa should be around 4.7-4.8
        assert 4.5 < pks["pka"] < 5.0
        # Ammonia pKn should be around 9.2-9.5
        assert 9.0 < pks["pkn"] < 9.8
        # Phosphate pKp should be around 7.0-7.3
        assert 6.8 < pks["pkp"] < 7.6

    def test_activity_corrected_pk_differ_from_uncorrected(self):
        from titration.chemistry import calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        # Activity-corrected values should differ from uncorrected
        assert pks["pk11"] != pks["pk1"]
        assert pks["pk22"] != pks["pk2"]
        assert pks["pkaa"] != pks["pka"]
        assert pks["pknn"] != pks["pkn"]
        assert pks["pkpp"] != pks["pkp"]

    def test_logf2_is_four_times_logf1(self):
        from titration.chemistry import calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        np.testing.assert_allclose(pks["logf2"], 4.0 * pks["logf1"], rtol=1e-10)

    def test_low_tds_adjustment(self):
        from titration.chemistry import calculate_pk_constants

        # TDS < 20 should be adjusted to 21 internally
        pks_low = calculate_pk_constants(21.0, 10.0, 5.0)
        pks_21 = calculate_pk_constants(21.0, 21.0, 5.0)
        np.testing.assert_allclose(pks_low["pk1"], pks_21["pk1"], rtol=1e-10)
        np.testing.assert_allclose(pks_low["pk11"], pks_21["pk11"], rtol=1e-10)

    def test_temperature_affects_pk(self):
        from titration.chemistry import calculate_pk_constants

        pks_low = calculate_pk_constants(10.0, 3300.0, 5.0)
        pks_high = calculate_pk_constants(35.0, 3300.0, 5.0)
        # pK values are temperature-dependent, so they should differ
        assert pks_low["pk1"] != pks_high["pk1"]
        assert pks_low["pk2"] != pks_high["pk2"]
        assert pks_low["pka"] != pks_high["pka"]


class TestSpeciesFractions:
    """Verify carbonate species fractions sum to 1.0."""

    @pytest.mark.parametrize("ph", [3.0, 4.0, 5.0, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0])
    def test_species_sum_to_one(self, ph):
        from titration.chemistry import per_h2co3, per_hco3, per_co3, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        pk11, pk22 = pks["pk11"], pks["pk22"]
        total = per_h2co3(ph, pk11, pk22) + per_hco3(ph, pk11, pk22) + per_co3(ph, pk11, pk22)
        np.testing.assert_allclose(total, 1.0, atol=1e-12)

    def test_h2co3_dominates_at_low_ph(self):
        from titration.chemistry import per_h2co3, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        frac = per_h2co3(3.0, pks["pk11"], pks["pk22"])
        assert frac > 0.99

    def test_hco3_dominates_at_mid_ph(self):
        from titration.chemistry import per_hco3, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        # At pH ~8.3 (between pK1 and pK2), HCO3- dominates
        frac = per_hco3(8.3, pks["pk11"], pks["pk22"])
        assert frac > 0.9

    def test_co3_dominates_at_high_ph(self):
        from titration.chemistry import per_co3, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        frac = per_co3(13.0, pks["pk11"], pks["pk22"])
        assert frac > 0.99

    def test_all_fractions_non_negative(self):
        from titration.chemistry import per_h2co3, per_hco3, per_co3, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        pk11, pk22 = pks["pk11"], pks["pk22"]
        for ph in [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]:
            assert per_h2co3(ph, pk11, pk22) >= 0.0
            assert per_hco3(ph, pk11, pk22) >= 0.0
            assert per_co3(ph, pk11, pk22) >= 0.0


class TestGenericWeakAcidFraction:
    """Verify the generic weak acid deprotonated fraction function."""

    def test_at_pka_equals_half(self):
        from titration.chemistry import per

        # At pH = pK, fraction should be exactly 0.5
        np.testing.assert_allclose(per(4.75, 4.75), 0.5, atol=1e-12)

    def test_below_pka_mostly_protonated(self):
        from titration.chemistry import per

        # pH well below pK: mostly protonated (fraction << 0.5)
        frac = per(2.0, 4.75)
        assert frac < 0.01

    def test_above_pka_mostly_deprotonated(self):
        from titration.chemistry import per

        # pH well above pK: mostly deprotonated (fraction >> 0.5)
        frac = per(8.0, 4.75)
        assert frac > 0.99

    def test_symmetric_around_pka(self):
        from titration.chemistry import per

        # per(pK+x) = 1 - per(pK-x) for symmetric Boltzmann distribution
        pka = 4.75
        delta = 1.5
        f_above = per(pka + delta, pka)
        f_below = per(pka - delta, pka)
        np.testing.assert_allclose(f_above + f_below, 1.0, atol=1e-12)


class TestAlkalinityDifferenceFunctions:
    """Verify d_h2co3_alk and d_hac_alk."""

    def test_d_h2co3_alk_same_ph_returns_zero(self):
        from titration.chemistry import d_h2co3_alk, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = d_h2co3_alk(6.0, 6.0, pks["pk11"], pks["pk22"])
        np.testing.assert_allclose(result, 0.0, atol=1e-15)

    def test_d_hac_alk_same_ph_returns_zero(self):
        from titration.chemistry import d_hac_alk, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = d_hac_alk(6.0, 6.0, pks["pkaa"])
        np.testing.assert_allclose(result, 0.0, atol=1e-15)

    def test_d_h2co3_alk_nonzero_for_different_ph(self):
        from titration.chemistry import d_h2co3_alk, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = d_h2co3_alk(6.75, 5.95, pks["pk11"], pks["pk22"])
        assert result != 0.0

    def test_d_hac_alk_nonzero_for_different_ph(self):
        from titration.chemistry import d_hac_alk, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = d_hac_alk(6.75, 5.95, pks["pkaa"])
        assert result != 0.0


class TestMassBalanceCorrections:
    """Verify m_h2o, m_nh3, m_hpo4 functions."""

    def test_m_h2o_nonzero_for_different_ph(self):
        from titration.chemistry import m_h2o, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = m_h2o(1.06, 3.50, 6.75, 5.95, 50.0, pks["logf1"])
        assert result != 0.0

    def test_m_nh3_zero_when_no_nitrogen(self):
        from titration.chemistry import m_nh3, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = m_nh3(6.75, 5.95, 0.0, 5.0, 50.0, pks["pknn"])
        np.testing.assert_allclose(result, 0.0, atol=1e-15)

    def test_m_nh3_nonzero_with_nitrogen(self):
        from titration.chemistry import m_nh3, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = m_nh3(6.75, 5.95, 100.0, 5.0, 50.0, pks["pknn"])
        assert result != 0.0

    def test_m_hpo4_zero_when_no_phosphorus(self):
        from titration.chemistry import m_hpo4, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = m_hpo4(6.75, 5.95, 0.0, 5.0, 50.0, pks["pkpp"])
        np.testing.assert_allclose(result, 0.0, atol=1e-15)

    def test_m_hpo4_nonzero_with_phosphorus(self):
        from titration.chemistry import m_hpo4, calculate_pk_constants

        pks = calculate_pk_constants(21.0, 3300.0, 5.0)
        result = m_hpo4(6.75, 5.95, 50.0, 5.0, 50.0, pks["pkpp"])
        assert result != 0.0
