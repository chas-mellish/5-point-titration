"""Functional tests for the titration CLI command."""

import subprocess
import sys
from pathlib import Path

WORKING_DIR = Path(__file__).resolve().parents[3]


def run_cli(*args, input_text=None):
    """Helper to invoke the CLI and capture output."""
    result = subprocess.run(
        [sys.executable, "-m", "titration.cli", *args],
        cwd=WORKING_DIR,
        capture_output=True,
        text=True,
        timeout=30,
        input=input_text,
    )
    return result


class TestCLIHelp:
    """CLI --help flag."""

    def test_help_output(self):
        result = run_cli("--help")
        assert result.returncode == 0
        assert "usage" in result.stdout.lower()
        assert "titration" in result.stdout.lower()

    def test_help_shows_all_parameters(self):
        result = run_cli("--help")
        assert result.returncode == 0
        assert "--ph0" in result.stdout
        assert "--ph1" in result.stdout
        assert "--ph2" in result.stdout
        assert "--ph3" in result.stdout
        assert "--ph4" in result.stdout
        assert "--vx1" in result.stdout
        assert "--vx2" in result.stdout
        assert "--vx3" in result.stdout
        assert "--vx4" in result.stdout
        assert "--titrant-normality" in result.stdout
        assert "--sample-volume-undiluted" in result.stdout
        assert "--sample-volume-diluted" in result.stdout
        assert "--temperature" in result.stdout
        assert "--tds" in result.stdout
        assert "--inorganic-nitrogen" in result.stdout
        assert "--inorganic-phosphorus" in result.stdout

    def test_help_shows_description(self):
        result = run_cli("--help")
        assert result.returncode == 0
        assert "5-point titration" in result.stdout.lower()


class TestCLIDefaults:
    """CLI with thesis default values (no arguments)."""

    def test_runs_with_no_args(self):
        result = run_cli()
        assert result.returncode == 0

    def test_output_contains_alkalinity(self):
        result = run_cli()
        assert result.returncode == 0
        assert "H2CO3* alkalinity" in result.stdout

    def test_output_contains_scfa(self):
        result = run_cli()
        assert result.returncode == 0
        assert "Short-chain fatty acids" in result.stdout

    def test_output_contains_ph_error(self):
        result = run_cli()
        assert result.returncode == 0
        assert "Systematic pH error" in result.stdout

    def test_output_contains_convergence(self):
        result = run_cli()
        assert result.returncode == 0
        assert "Convergence: OK" in result.stdout

    def test_output_contains_results_header(self):
        result = run_cli()
        assert result.returncode == 0
        assert "5-Point Titration Results" in result.stdout

    def test_output_alkalinity_value(self):
        result = run_cli()
        assert result.returncode == 0
        assert "1863.85" in result.stdout

    def test_output_scfa_value(self):
        result = run_cli()
        assert result.returncode == 0
        assert "195.91" in result.stdout


class TestCLICustomParameters:
    """CLI with custom parameter values."""

    def test_custom_temperature(self):
        result = run_cli("--temperature", "35.0")
        assert result.returncode == 0
        assert "H2CO3* alkalinity" in result.stdout
        assert "Convergence: OK" in result.stdout

    def test_custom_ph_values(self):
        result = run_cli("--ph0", "7.5", "--ph1", "6.8", "--ph2", "6.0",
                         "--ph3", "5.2", "--ph4", "4.3")
        assert result.returncode == 0
        assert "H2CO3* alkalinity" in result.stdout

    def test_custom_tds(self):
        result = run_cli("--tds", "5000.0")
        assert result.returncode == 0
        assert "H2CO3* alkalinity" in result.stdout

    def test_with_nitrogen(self):
        result = run_cli("--inorganic-nitrogen", "100.0")
        assert result.returncode == 0
        assert "H2CO3* alkalinity" in result.stdout

    def test_with_phosphorus(self):
        result = run_cli("--inorganic-phosphorus", "50.0")
        assert result.returncode == 0
        assert "H2CO3* alkalinity" in result.stdout


class TestCLIInvalidInput:
    """CLI with invalid arguments."""

    def test_unknown_flag(self):
        result = run_cli("--nonexistent-flag")
        assert result.returncode != 0

    def test_non_numeric_ph(self):
        result = run_cli("--ph0", "not_a_number")
        assert result.returncode != 0

    def test_non_numeric_temperature(self):
        result = run_cli("--temperature", "abc")
        assert result.returncode != 0
