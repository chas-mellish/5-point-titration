import pytest

from titration.chemistry import calculate_pk_constants


@pytest.fixture
def thesis_pk_constants():
    return calculate_pk_constants(21.0, 3300.0, 5.0)
