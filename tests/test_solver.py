"""Integration tests for titration.solver — run_solver pipeline."""

import dataclasses

import numpy as np
import pytest

from titration.constants import THESIS_DEFAULTS
from titration.models import TitrationInput
from titration.solver import run_solver


class TestThesisDefaults:
    """Verify solver output with the canonical thesis default inputs."""

    def test_convergence_status(self):
        """Thesis defaults must converge."""
        result = run_solver(THESIS_DEFAULTS)
        assert result.convergence_status == "converged"

    def test_h2co3_alkalinity_positive_and_reasonable(self):
        """H2CO3* alkalinity should be positive and in the hundreds-to-thousands
        mg/L range for typical anaerobic-digester effluent."""
        result = run_solver(THESIS_DEFAULTS)
        assert result.h2co3_alkalinity > 0
        assert 100 < result.h2co3_alkalinity < 10000

    def test_scfa_concentration_non_negative(self):
        """SCFA concentration is clamped to >= 0 by the solver."""
        result = run_solver(THESIS_DEFAULTS)
        assert result.scfa_concentration >= 0.0

    def test_systematic_ph_error_small(self):
        """Systematic pH error should be small (|error| < 0.2)."""
        result = run_solver(THESIS_DEFAULTS)
        assert abs(result.systematic_ph_error) < 0.2

    def test_ct_values_close_when_converged(self):
        """When converged, the two independent Ct estimates should be close."""
        result = run_solver(THESIS_DEFAULTS)
        ct1, ct2 = result.ct_values
        np.testing.assert_allclose(ct1, ct2, rtol=0.01)

    def test_at_values_are_tuple(self):
        """at_values should be a two-element tuple."""
        result = run_solver(THESIS_DEFAULTS)
        assert len(result.at_values) == 2


class TestZeroNitrogenPhosphorus:
    """Verify behaviour when Nt=0 and Pt=0 (carbonate-only system)."""

    def test_zero_nt_pt_converges(self):
        """With Nt=0 and Pt=0, solver still converges (thesis defaults
        already have these zeroed)."""
        result = run_solver(THESIS_DEFAULTS)
        assert result.convergence_status == "converged"

    def test_results_unchanged_with_explicit_zero_nutrients(self):
        """Explicit zero Nt/Pt should match the default (which is also zero)."""
        inp = dataclasses.replace(
            THESIS_DEFAULTS, inorganic_nitrogen=0.0, inorganic_phosphorus=0.0,
        )
        r1 = run_solver(THESIS_DEFAULTS)
        r2 = run_solver(inp)
        np.testing.assert_allclose(
            r1.h2co3_alkalinity, r2.h2co3_alkalinity, rtol=1e-12,
        )


class TestEdgeCaseLowTDS:
    """TDS below the 20-threshold should be clamped to 21."""

    def test_low_tds_undiluted_converges(self):
        """TDS=15 with dil=1 (undiluted) should converge without error."""
        inp = dataclasses.replace(
            THESIS_DEFAULTS, sample_volume_undiluted=50.0, tds=15.0,
        )
        result = run_solver(inp)
        assert result.convergence_status == "converged"
        assert np.isfinite(result.h2co3_alkalinity)

    def test_low_tds_matches_tds21(self):
        """TDS=15 should be clamped to 21, so results must match TDS=21."""
        r15 = run_solver(dataclasses.replace(THESIS_DEFAULTS, sample_volume_undiluted=50.0, tds=15.0))
        r21 = run_solver(dataclasses.replace(THESIS_DEFAULTS, sample_volume_undiluted=50.0, tds=21.0))
        np.testing.assert_allclose(
            r15.h2co3_alkalinity, r21.h2co3_alkalinity, rtol=1e-12,
        )


class TestExtremeTemperature:
    """Solver must handle cold and warm temperatures gracefully."""

    def test_cold_temperature_converges(self):
        """Temperature=5 C (cold water) should converge."""
        inp = dataclasses.replace(THESIS_DEFAULTS, temperature=5.0)
        result = run_solver(inp)
        assert result.convergence_status == "converged"
        assert result.h2co3_alkalinity > 0

    def test_warm_temperature_converges(self):
        """Temperature=40 C (warm water) should converge."""
        inp = dataclasses.replace(THESIS_DEFAULTS, temperature=40.0)
        result = run_solver(inp)
        assert result.convergence_status == "converged"
        assert result.h2co3_alkalinity > 0

    def test_temperature_affects_results(self):
        """Different temperatures should produce different alkalinity values."""
        r_cold = run_solver(dataclasses.replace(THESIS_DEFAULTS, temperature=5.0))
        r_warm = run_solver(dataclasses.replace(THESIS_DEFAULTS, temperature=40.0))
        assert r_cold.h2co3_alkalinity != r_warm.h2co3_alkalinity


class TestAdversarialInputs:
    """Verify division-by-zero guards with adversarial pH values."""

    _base = dict(
        vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
        titrant_normality=0.0728,
        sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
        temperature=21.0, tds=3300.0,
    )

    def test_identical_ph3_ph4_raises(self):
        """ph3 == ph4 causes d_hac_alk denominator to be zero."""
        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=5.18,
            **self._base,
        )
        with pytest.raises(ValueError, match="ph3 and ph4 must differ"):
            run_solver(inp)

    def test_identical_ph1_ph4_raises(self):
        """ph1 == ph4 causes d_hac_alk(ph1, ph4) denominator to be zero."""
        inp = TitrationInput(
            ph0=7.36, ph1=4.29, ph2=5.95, ph3=5.18, ph4=4.29,
            **self._base,
        )
        with pytest.raises(ValueError, match="ph1 and ph4 must differ"):
            run_solver(inp)


class TestInputValidation:
    """Verify sample volume validation guards."""

    def test_zero_undiluted_raises(self):
        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=0.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
        )
        with pytest.raises(ValueError, match="sample_volume_undiluted"):
            run_solver(inp)

    def test_negative_undiluted_raises(self):
        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=-10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
        )
        with pytest.raises(ValueError, match="sample_volume_undiluted"):
            run_solver(inp)

    def test_zero_diluted_raises(self):
        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=0.0,
            temperature=21.0, tds=3300.0,
        )
        with pytest.raises(ValueError, match="sample_volume_diluted"):
            run_solver(inp)
