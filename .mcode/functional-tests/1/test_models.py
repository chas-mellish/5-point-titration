"""Functional tests for titration.models — TitrationInput and TitrationResult dataclasses."""

import pytest


class TestTitrationInputCreation:
    """Verify TitrationInput dataclass construction with all 16 parameters."""

    def test_create_with_all_required_fields(self):
        from titration.models import TitrationInput

        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
        )
        assert inp.ph0 == 7.36
        assert inp.ph1 == 6.75
        assert inp.ph2 == 5.95
        assert inp.ph3 == 5.18
        assert inp.ph4 == 4.29
        assert inp.vx1 == 1.06
        assert inp.vx2 == 3.50
        assert inp.vx3 == 4.84
        assert inp.vx4 == 5.40
        assert inp.titrant_normality == 0.0728
        assert inp.sample_volume_undiluted == 10.0
        assert inp.sample_volume_diluted == 50.0
        assert inp.temperature == 21.0
        assert inp.tds == 3300.0
        assert inp.inorganic_nitrogen == 0.0
        assert inp.inorganic_phosphorus == 0.0

    def test_optional_nitrogen_defaults_to_zero(self):
        from titration.models import TitrationInput

        inp = TitrationInput(
            ph0=7.0, ph1=6.5, ph2=5.5, ph3=5.0, ph4=4.0,
            vx1=1.0, vx2=3.0, vx3=4.5, vx4=5.0,
            titrant_normality=0.1,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=25.0, tds=1000.0,
        )
        assert inp.inorganic_nitrogen == 0.0

    def test_optional_phosphorus_defaults_to_zero(self):
        from titration.models import TitrationInput

        inp = TitrationInput(
            ph0=7.0, ph1=6.5, ph2=5.5, ph3=5.0, ph4=4.0,
            vx1=1.0, vx2=3.0, vx3=4.5, vx4=5.0,
            titrant_normality=0.1,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=25.0, tds=1000.0,
        )
        assert inp.inorganic_phosphorus == 0.0

    def test_create_with_nitrogen_and_phosphorus(self):
        from titration.models import TitrationInput

        inp = TitrationInput(
            ph0=7.36, ph1=6.75, ph2=5.95, ph3=5.18, ph4=4.29,
            vx1=1.06, vx2=3.50, vx3=4.84, vx4=5.40,
            titrant_normality=0.0728,
            sample_volume_undiluted=10.0, sample_volume_diluted=50.0,
            temperature=21.0, tds=3300.0,
            inorganic_nitrogen=100.0, inorganic_phosphorus=50.0,
        )
        assert inp.inorganic_nitrogen == 100.0
        assert inp.inorganic_phosphorus == 50.0

    def test_missing_required_field_raises_error(self):
        from titration.models import TitrationInput

        with pytest.raises(TypeError):
            TitrationInput(
                ph0=7.36, ph1=6.75,
                # Missing most required fields
            )


class TestTitrationInputFromDict:
    """Verify TitrationInput.from_dict() classmethod."""

    def test_from_dict_with_all_fields(self):
        from titration.models import TitrationInput

        data = {
            "ph0": 7.36, "ph1": 6.75, "ph2": 5.95, "ph3": 5.18, "ph4": 4.29,
            "vx1": 1.06, "vx2": 3.50, "vx3": 4.84, "vx4": 5.40,
            "titrant_normality": 0.0728,
            "sample_volume_undiluted": 10.0, "sample_volume_diluted": 50.0,
            "temperature": 21.0, "tds": 3300.0,
            "inorganic_nitrogen": 100.0, "inorganic_phosphorus": 50.0,
        }
        inp = TitrationInput.from_dict(data)
        assert inp.ph0 == 7.36
        assert inp.temperature == 21.0
        assert inp.inorganic_nitrogen == 100.0
        assert inp.inorganic_phosphorus == 50.0

    def test_from_dict_ignores_extra_fields(self):
        from titration.models import TitrationInput

        data = {
            "ph0": 7.36, "ph1": 6.75, "ph2": 5.95, "ph3": 5.18, "ph4": 4.29,
            "vx1": 1.06, "vx2": 3.50, "vx3": 4.84, "vx4": 5.40,
            "titrant_normality": 0.0728,
            "sample_volume_undiluted": 10.0, "sample_volume_diluted": 50.0,
            "temperature": 21.0, "tds": 3300.0,
            "extra_field": "should be ignored",
            "another_extra": 42,
        }
        inp = TitrationInput.from_dict(data)
        assert inp.ph0 == 7.36
        assert not hasattr(inp, "extra_field")

    def test_from_dict_with_optional_defaults(self):
        from titration.models import TitrationInput

        data = {
            "ph0": 7.36, "ph1": 6.75, "ph2": 5.95, "ph3": 5.18, "ph4": 4.29,
            "vx1": 1.06, "vx2": 3.50, "vx3": 4.84, "vx4": 5.40,
            "titrant_normality": 0.0728,
            "sample_volume_undiluted": 10.0, "sample_volume_diluted": 50.0,
            "temperature": 21.0, "tds": 3300.0,
        }
        inp = TitrationInput.from_dict(data)
        assert inp.inorganic_nitrogen == 0.0
        assert inp.inorganic_phosphorus == 0.0

    def test_from_dict_missing_required_raises_error(self):
        from titration.models import TitrationInput

        data = {"ph0": 7.36, "ph1": 6.75}
        with pytest.raises(TypeError):
            TitrationInput.from_dict(data)


class TestTitrationResult:
    """Verify TitrationResult dataclass structure."""

    def test_create_result(self):
        from titration.models import TitrationResult

        result = TitrationResult(
            h2co3_alkalinity=1863.85,
            scfa_concentration=195.91,
            systematic_ph_error=-0.03,
            convergence_status="converged",
            ct_values=(2043.04, 2043.23),
            at_values=(195.91, 196.12),
        )
        assert result.h2co3_alkalinity == 1863.85
        assert result.scfa_concentration == 195.91
        assert result.systematic_ph_error == -0.03
        assert result.convergence_status == "converged"
        assert len(result.ct_values) == 2
        assert len(result.at_values) == 2
