import pytest
import numpy as np

@pytest.fixture
def sample_ph_data():
    """Fixture: Sample pH titration curve data"""
    return np.array([2.1, 2.5, 3.2, 4.1, 5.8, 7.2, 8.9, 9.5])

@pytest.fixture
def sample_volume_data():
    """Fixture: Corresponding titrant volumes (mL)"""
    return np.array([0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5])

@pytest.fixture
def expected_alkalinity():
    """Fixture: Expected alkalinity result (mg/L as CaCO3)"""
    return 125.5  # Example value from thesis