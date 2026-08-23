from __future__ import annotations

from .models import TitrationInput

CACO3_FACTOR = 50000
ACETIC_ACID_FACTOR = 60000
NITROGEN_FACTOR = 14000
PHOSPHORUS_FACTOR = 31000

THESIS_DEFAULTS = TitrationInput(
    ph0=7.36,
    ph1=6.75,
    ph2=5.95,
    ph3=5.18,
    ph4=4.29,
    vx1=1.06,
    vx2=3.50,
    vx3=4.84,
    vx4=5.40,
    titrant_normality=0.0728,
    sample_volume_undiluted=10.0,
    sample_volume_diluted=50.0,
    temperature=21.0,
    tds=3300.0,
    inorganic_nitrogen=0.0,
    inorganic_phosphorus=0.0,
)
