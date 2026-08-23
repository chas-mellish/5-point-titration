from __future__ import annotations

from dataclasses import dataclass


@dataclass
class TitrationInput:
    ph0: float
    ph1: float
    ph2: float
    ph3: float
    ph4: float
    vx1: float
    vx2: float
    vx3: float
    vx4: float
    titrant_normality: float
    sample_volume_undiluted: float
    sample_volume_diluted: float
    temperature: float
    tds: float
    inorganic_nitrogen: float = 0.0
    inorganic_phosphorus: float = 0.0

    @classmethod
    def from_dict(cls, data: dict) -> TitrationInput:
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class TitrationResult:
    h2co3_alkalinity: float
    scfa_concentration: float
    systematic_ph_error: float
    convergence_status: str
    ct_values: tuple[float, float]
    at_values: tuple[float, float]
