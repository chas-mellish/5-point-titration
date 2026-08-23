"""Functional tests for titration.solver and titration.core — the solver pipeline and orchestrator."""

import numpy as np


class TestRunTitrationThesisDefaults:
    """Verify run_titration with thesis default values produces converged results."""

    def test_thesis_defaults_converges(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        assert result.convergence_status == "converged"

    def test_thesis_defaults_alkalinity_positive(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        assert result.h2co3_alkalinity > 0.0

    def test_thesis_defaults_alkalinity_reasonable_range(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        # H2CO3* alkalinity for anaerobic digester should be in the hundreds to
        # low thousands of mg/L as CaCO3
        assert 500.0 < result.h2co3_alkalinity < 5000.0

    def test_thesis_defaults_alkalinity_value(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        np.testing.assert_allclose(result.h2co3_alkalinity, 1863.85, rtol=0.01)

    def test_thesis_defaults_scfa_non_negative(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        assert result.scfa_concentration >= 0.0

    def test_thesis_defaults_scfa_value(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        np.testing.assert_allclose(result.scfa_concentration, 195.91, rtol=0.01)

    def test_thesis_defaults_ph_error_small(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        # Systematic pH error should be small (within +/- 0.2)
        assert abs(result.systematic_ph_error) < 0.2

    def test_thesis_defaults_ph_error_value(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        np.testing.assert_allclose(result.systematic_ph_error, -0.03, atol=0.005)

    def test_thesis_defaults_ct_values_close(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        ct1, ct2 = result.ct_values
        # When converged, Ct1 and Ct2 should be close (sign change means they bracket)
        assert abs(ct1 - ct2) < 10.0

    def test_thesis_defaults_at_values_tuple(self):
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        assert len(result.at_values) == 2
        assert isinstance(result.at_values[0], (float, np.floating))
        assert isinstance(result.at_values[1], (float, np.floating))


class TestRunTitrationWithNitrogenPhosphorus:
    """Verify run_titration with nitrogen and phosphorus corrections."""

    def test_with_nitrogen_converges(self):
        from titration import run_titration, TitrationInput

        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
            inorganic_nitrogen=100.0, inorganic_phosphorus=0.0,
        )
        result = run_titration(inp)
        assert result.convergence_status == "converged"

    def test_with_phosphorus_converges(self):
        from titration import run_titration, TitrationInput

        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
            inorganic_nitrogen=0.0, inorganic_phosphorus=50.0,
        )
        result = run_titration(inp)
        assert result.convergence_status == "converged"

    def test_with_both_np_converges(self):
        from titration import run_titration, TitrationInput

        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
            inorganic_nitrogen=100.0, inorganic_phosphorus=50.0,
        )
        result = run_titration(inp)
        assert result.convergence_status == "converged"
        np.testing.assert_allclose(result.h2co3_alkalinity, 1822.90, rtol=0.01)

    def test_nitrogen_affects_alkalinity(self):
        from titration import run_titration, TitrationInput, THESIS_DEFAULTS

        result_no_n = run_titration(THESIS_DEFAULTS)

        inp_with_n = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
            inorganic_nitrogen=100.0,
        )
        result_with_n = run_titration(inp_with_n)
        assert result_no_n.h2co3_alkalinity != result_with_n.h2co3_alkalinity


class TestRunTitrationTemperatureVariation:
    """Verify run_titration at different temperatures."""

    def test_higher_temperature_converges(self):
        from titration import run_titration, TitrationInput

        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=35.0, tds=3300.0,
        )
        result = run_titration(inp)
        assert result.convergence_status == "converged"
        assert result.h2co3_alkalinity > 0.0

    def test_lower_temperature_converges(self):
        from titration import run_titration, TitrationInput

        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=10.0, tds=3300.0,
        )
        result = run_titration(inp)
        assert result.convergence_status == "converged"
        assert result.h2co3_alkalinity > 0.0

    def test_temperature_affects_results(self):
        from titration import run_titration, TitrationInput

        def make_input(temp):
            return TitrationInput(
                ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
                vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
                titrant_normality=0.0728,
                sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
                temperature=temp, tds=3300.0,
            )

        r_10 = run_titration(make_input(10.0))
        r_35 = run_titration(make_input(35.0))
        # Different temperatures should produce different alkalinity values
        assert r_10.h2co3_alkalinity != r_35.h2co3_alkalinity


class TestSolverEdgeCases:
    """Verify solver convergence edge cases."""

    def test_close_ph_values_exceed_iterations(self):
        from titration import run_titration, TitrationInput

        # Very close pH values should cause convergence issues
        inp = TitrationInput(
            ph0=7.00, ph1=6.99, ph2=6.98, ph3=6.97, ph4=6.96,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
        )
        result = run_titration(inp)
        # Should report a non-converged status (exceeded iterations or ratio too high)
        assert result.convergence_status in ("exceeded_max_iterations", "ratio_too_high")

    def test_high_tds_converges(self):
        from titration import run_titration, TitrationInput

        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=10000.0,
        )
        result = run_titration(inp)
        assert result.convergence_status == "converged"
        assert result.h2co3_alkalinity > 0.0

    def test_run_titration_delegates_to_solver(self):
        """run_titration in core.py should delegate to run_solver and return same result."""
        from titration.core import run_titration
        from titration.solver import run_solver
        from titration import THESIS_DEFAULTS

        r1 = run_titration(THESIS_DEFAULTS)
        r2 = run_solver(THESIS_DEFAULTS)
        np.testing.assert_allclose(r1.h2co3_alkalinity, r2.h2co3_alkalinity, rtol=1e-10)
        np.testing.assert_allclose(r1.scfa_concentration, r2.scfa_concentration, rtol=1e-10)
        assert r1.convergence_status == r2.convergence_status

    def test_scfa_clamped_to_zero(self):
        """SCFA should never be negative (clamped to 0 per Pascal output procedure)."""
        from titration import run_titration, THESIS_DEFAULTS

        result = run_titration(THESIS_DEFAULTS)
        assert result.scfa_concentration >= 0.0


class TestConstants:
    """Verify conversion factor constants match spec."""

    def test_caco3_factor(self):
        from titration.constants import CACO3_FACTOR
        assert CACO3_FACTOR == 50000

    def test_acetic_acid_factor(self):
        from titration.constants import ACETIC_ACID_FACTOR
        assert ACETIC_ACID_FACTOR == 60000

    def test_nitrogen_factor(self):
        from titration.constants import NITROGEN_FACTOR
        assert NITROGEN_FACTOR == 14000

    def test_phosphorus_factor(self):
        from titration.constants import PHOSPHORUS_FACTOR
        assert PHOSPHORUS_FACTOR == 31000

    def test_thesis_defaults_values(self):
        from titration.constants import THESIS_DEFAULTS

        assert THESIS_DEFAULTS.ph0 == 7.36
        assert THESIS_DEFAULTS.ph1 == 6.75
        assert THESIS_DEFAULTS.ph2 == 5.95
        assert THESIS_DEFAULTS.ph3 == 5.18
        assert THESIS_DEFAULTS.ph4 == 4.29
        assert THESIS_DEFAULTS.vx1 == 1.06
        assert THESIS_DEFAULTS.vx2 == 3.50
        assert THESIS_DEFAULTS.vx3 == 4.84
        assert THESIS_DEFAULTS.vx4 == 5.40
        assert THESIS_DEFAULTS.titrant_normality == 0.0728
        assert THESIS_DEFAULTS.sample_volume_undiluted == 10.0
        assert THESIS_DEFAULTS.sample_volume_diluted == 50.0
        assert THESIS_DEFAULTS.temperature == 21.0
        assert THESIS_DEFAULTS.tds == 3300.0
        assert THESIS_DEFAULTS.inorganic_nitrogen == 0.0
        assert THESIS_DEFAULTS.inorganic_phosphorus == 0.0


class TestPublicApiExports:
    """Verify __init__.py exports match spec."""

    def test_run_titration_importable(self):
        from titration import run_titration
        assert callable(run_titration)

    def test_titration_input_importable(self):
        from titration import TitrationInput
        assert TitrationInput is not None

    def test_titration_result_importable(self):
        from titration import TitrationResult
        assert TitrationResult is not None

    def test_thesis_defaults_importable(self):
        from titration import THESIS_DEFAULTS
        assert THESIS_DEFAULTS is not None

    def test_version_importable(self):
        from titration import __version__
        assert isinstance(__version__, str)
        assert __version__ == "0.1.0"
