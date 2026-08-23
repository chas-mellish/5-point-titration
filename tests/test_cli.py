"""Tests for the titration CLI module."""

import argparse

import pytest

from titration.cli import _build_parser, _format_results, main
from titration.constants import THESIS_DEFAULTS
from titration.models import TitrationResult


class TestBuildParser:
    """Verify argparse configuration."""

    def test_returns_argument_parser(self):
        p = _build_parser()
        assert isinstance(p, argparse.ArgumentParser)

    def test_defaults_match_thesis_defaults(self):
        p = _build_parser()
        args = p.parse_args([])
        assert args.ph0 == THESIS_DEFAULTS.ph0
        assert args.ph4 == THESIS_DEFAULTS.ph4
        assert args.titrant_normality == THESIS_DEFAULTS.titrant_normality
        assert args.temperature == THESIS_DEFAULTS.temperature
        assert args.tds == THESIS_DEFAULTS.tds

    def test_custom_values_parse(self):
        p = _build_parser()
        args = p.parse_args(["--ph0", "8.0", "--temperature", "30.0"])
        assert args.ph0 == 8.0
        assert args.temperature == 30.0


class TestFormatResults:
    """Verify output formatting."""

    def _make_result(self, status="converged"):
        return TitrationResult(
            h2co3_alkalinity=2500.0,
            scfa_concentration=150.0,
            systematic_ph_error=-0.02,
            convergence_status=status,
            ct_values=(100.0, 99.5),
            at_values=(150.0, 145.0),
        )

    def test_converged_output_contains_ok(self):
        text = _format_results(self._make_result("converged"))
        assert "Convergence: OK" in text

    def test_exceeded_iterations_warning(self):
        text = _format_results(self._make_result("exceeded_max_iterations"))
        assert "Maximum iterations exceeded" in text

    def test_ratio_too_high_warning(self):
        text = _format_results(self._make_result("ratio_too_high"))
        assert "ratio too high" in text

    def test_alkalinity_in_output(self):
        text = _format_results(self._make_result())
        assert "2500.00" in text

    def test_negative_scfa_clamped(self):
        r = TitrationResult(
            h2co3_alkalinity=2500.0,
            scfa_concentration=-50.0,
            systematic_ph_error=0.0,
            convergence_status="converged",
            ct_values=(100.0, 100.0),
            at_values=(-50.0, -45.0),
        )
        text = _format_results(r)
        assert "-50" not in text


class TestMain:
    """Verify the main entry point."""

    def test_main_with_defaults(self, capsys):
        main([])
        captured = capsys.readouterr()
        assert "H2CO3* alkalinity" in captured.out
        assert "Convergence: OK" in captured.out

    def test_main_with_custom_args(self, capsys):
        main(["--temperature", "30.0", "--tds", "5000.0"])
        captured = capsys.readouterr()
        assert "H2CO3* alkalinity" in captured.out

    def test_main_invalid_arg_exits(self):
        with pytest.raises(SystemExit):
            main(["--invalid-flag"])
