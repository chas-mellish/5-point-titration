"""Parametric and boundary tests for the titration public API."""

import dataclasses

import numpy as np
import pytest

from titration.constants import THESIS_DEFAULTS
from titration.core import run_titration
from titration.models import TitrationInput, TitrationResult
from titration.solver import run_solver


# ---------------------------------------------------------------------------
# Parametric tests
# ---------------------------------------------------------------------------

class TestParametricTemperature:
    """Run the solver across a range of temperatures."""

    @pytest.mark.parametrize("temperature", [5, 10, 15, 20, 25, 30, 35, 40])
    def test_solver_converges_at_temperature(self, temperature):
        """Solver should converge for each reasonable temperature value."""
        inp = dataclasses.replace(THESIS_DEFAULTS, temperature=float(temperature))
        result = run_solver(inp)
        assert result.convergence_status == "converged"
        assert result.h2co3_alkalinity > 0


class TestParametricTDS:
    """Run the solver across a range of TDS values."""

    @pytest.mark.parametrize("tds", [500, 1000, 2000, 3300, 5000, 8000])
    def test_solver_converges_at_tds(self, tds):
        """Solver should converge for each TDS value."""
        inp = dataclasses.replace(THESIS_DEFAULTS, tds=float(tds))
        result = run_solver(inp)
        assert result.convergence_status == "converged"
        assert result.h2co3_alkalinity > 0


# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------

class TestBoundaryConditions:
    """Edge-case and boundary inputs."""

    def test_very_high_tds(self):
        """TDS=10000 should still converge and produce finite results."""
        inp = dataclasses.replace(THESIS_DEFAULTS, tds=10000.0)
        result = run_solver(inp)
        assert result.convergence_status == "converged"
        assert np.isfinite(result.h2co3_alkalinity)

    def test_tds_near_threshold_clamped(self):
        """TDS=19 (< 20) should be clamped to 21 and match TDS=21."""
        r19 = run_solver(dataclasses.replace(THESIS_DEFAULTS, sample_volume_undiluted=50.0, tds=19.0))
        r21 = run_solver(dataclasses.replace(THESIS_DEFAULTS, sample_volume_undiluted=50.0, tds=21.0))
        np.testing.assert_allclose(
            r19.h2co3_alkalinity, r21.h2co3_alkalinity, rtol=1e-12,
        )


# ---------------------------------------------------------------------------
# Convergence guards
# ---------------------------------------------------------------------------

class TestConvergenceGuards:
    """Verify that the solver correctly reports non-convergence."""

    def test_exceeded_max_iterations(self):
        """Tiny titrant volumes can produce ct_comp that never crosses zero,
        exhausting the 20-iteration limit."""
        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=0.01, vx2=0.02, vx3=0.03, vx4=0.04,
            titrant_normality=0.0728,
            sample_volume_undiluted=50.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
        )
        result = run_solver(inp)
        assert result.convergence_status == "exceeded_max_iterations"

    def test_ratio_too_high(self):
        """When SCFA-to-carbonate ratio > 0.5, solver reports ratio_too_high."""
        inp = TitrationInput(
            ph0=5.50, ph1=5.40, ph2=5.30, ph3=5.20, ph4=5.10,
            vx1=0.5, vx2=1.0, vx3=1.5, vx4=2.0,
            titrant_normality=0.0728,
            sample_volume_undiluted=50.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
        )
        result = run_solver(inp)
        assert result.convergence_status == "ratio_too_high"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

class TestPublicAPI:
    """Tests for the public-facing run_titration and model helpers."""

    def test_run_titration_matches_run_solver(self):
        """run_titration() is a thin wrapper and must produce identical
        results to run_solver()."""
        r_solver = run_solver(THESIS_DEFAULTS)
        r_api = run_titration(THESIS_DEFAULTS)
        assert r_solver == r_api

    def test_from_dict_round_trip(self):
        """TitrationInput.from_dict(asdict(inp)) should reconstruct the
        original dataclass."""
        d = dataclasses.asdict(THESIS_DEFAULTS)
        reconstructed = TitrationInput.from_dict(d)
        assert reconstructed == THESIS_DEFAULTS

    def test_from_dict_ignores_extra_keys(self):
        """from_dict should silently drop keys not in the dataclass."""
        d = dataclasses.asdict(THESIS_DEFAULTS)
        d["extra_key"] = 42
        reconstructed = TitrationInput.from_dict(d)
        assert reconstructed == THESIS_DEFAULTS

    def test_titration_result_fields(self):
        """TitrationResult should expose all documented attributes."""
        result = run_solver(THESIS_DEFAULTS)
        assert isinstance(result, TitrationResult)
        assert hasattr(result, "h2co3_alkalinity")
        assert hasattr(result, "scfa_concentration")
        assert hasattr(result, "systematic_ph_error")
        assert hasattr(result, "convergence_status")
        assert hasattr(result, "ct_values")
        assert hasattr(result, "at_values")
