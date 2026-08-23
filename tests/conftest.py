"""Shared fixtures for the 5-point titration test suite."""

import pytest

from titration.chemistry import calculate_pk_constants
from titration.constants import THESIS_DEFAULTS


@pytest.fixture
def thesis_defaults():
    """Return the canonical thesis default TitrationInput."""
    return THESIS_DEFAULTS


@pytest.fixture
def thesis_pk_constants():
    """Return pK constants computed at thesis conditions (21 C, TDS=3300, dil=5)."""
    return calculate_pk_constants(21.0, 3300.0, 5.0)
