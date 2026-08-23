from __future__ import annotations

import numpy as np


def calculate_ionic_strength(tds: float, dil: float) -> float:
    # Pascal: nue(TDS,dil) -- mue := 0.000025 * (TDS/dil - 20)
    return 0.000025 * (tds / dil - 20.0)


def calculate_log_activity(mue: float, ktemp: float) -> float:
    # Extended Debye-Huckel activity coefficient for monovalent ions
    return (
        -1.825e6
        * (78.3 * ktemp) ** (-1.5)
        * (np.sqrt(mue) / (1.0 + np.sqrt(mue)) - 0.3 * mue)
    )


def calculate_pk_constants(temperature: float, tds: float, dil: float) -> dict:
    ktemp = 273.0 + temperature
    tds_adj = tds if tds >= 20.0 else 21.0

    mue = calculate_ionic_strength(tds_adj, dil)
    logf1 = calculate_log_activity(mue, ktemp)
    logf2 = 4.0 * logf1

    # Plummer & Busenberg carbonate equilibrium constants
    pk1 = -1.0 * (
        -356.3094
        - 0.06091964 * ktemp
        + 21834.37 / ktemp
        + 126.8339 * np.log10(ktemp)
        - 1684915.0 / (ktemp * ktemp)
    )
    pk11 = pk1 + logf1

    pk2 = -1.0 * (
        -107.8871
        - 0.03252849 * ktemp
        + 5151.79 / ktemp
        + 38.92561 * np.log10(ktemp)
        - 563713.9 / (ktemp * ktemp)
    )
    pk22 = pk2 - logf1 + logf2

    # Acetic acid dissociation
    pka = 1170.5 / ktemp - 3.165 + 0.0134 * ktemp
    pkaa = pka + logf1

    # Ammonium dissociation
    pkn = 2835.8 / ktemp - 0.6322 + 0.00123 * ktemp
    pknn = pkn + logf1

    # Phosphate second dissociation
    pkp = 1979.5 / ktemp - 5.3541 + 0.01984 * ktemp
    pkpp = pkp - logf1 + logf2

    return {
        "pk1": pk1,
        "pk2": pk2,
        "pka": pka,
        "pkn": pkn,
        "pkp": pkp,
        "pk11": pk11,
        "pk22": pk22,
        "pkaa": pkaa,
        "pknn": pknn,
        "pkpp": pkpp,
        "logf1": logf1,
        "logf2": logf2,
    }


def per_h2co3(ph: float, pk11: float, pk22: float) -> float:
    # Fraction of total carbonate as H2CO3*
    return 1.0 / (
        1.0
        + 10.0 ** (ph - pk11)
        + 10.0 ** (2.0 * ph - pk11 - pk22)
    )


def per_hco3(ph: float, pk11: float, pk22: float) -> float:
    # Fraction of total carbonate as HCO3-
    return 1.0 / (
        10.0 ** (pk11 - ph)
        + 1.0
        + 10.0 ** (ph - pk22)
    )


def per_co3(ph: float, pk11: float, pk22: float) -> float:
    # Fraction of total carbonate as CO3^2-
    return 1.0 / (
        10.0 ** (pk11 + pk22 - 2.0 * ph)
        + 10.0 ** (pk22 - ph)
        + 1.0
    )


def per(ph: float, pkk: float) -> float:
    # Generic weak acid deprotonated fraction: A- / (HA + A-)
    return 1.0 / (1.0 + 10.0 ** (pkk - ph))


def d_h2co3_alk(ph_f: float, ph_s: float, pk11: float, pk22: float) -> float:
    return (
        per_hco3(ph_f, pk11, pk22) - per_hco3(ph_s, pk11, pk22)
        + 2.0 * (per_co3(ph_f, pk11, pk22) - per_co3(ph_s, pk11, pk22))
    )


def d_hac_alk(ph_f: float, ph_s: float, pkaa: float) -> float:
    return per(ph_f, pkaa) - per(ph_s, pkaa)


def m_h2o(
    vxfi: float,
    vxs: float,
    ph_fi: float,
    ph_s: float,
    vsdil: float,
    logf1: float,
) -> float:
    # Water mass balance correction: H+ and OH- contributions
    return (
        (vsdil + vxs) * 10.0 ** (-ph_s) / 10.0 ** logf1
        - (vsdil + vxfi) * 10.0 ** (-ph_fi) / 10.0 ** logf1
        + (vsdil + vxfi) * 10.0 ** (ph_fi - 14.0)
        - (vsdil + vxs) * 10.0 ** (ph_s - 14.0)
    )


def m_nh3(
    ph_f: float,
    ph_s: float,
    nt: float,
    dil: float,
    vsdil: float,
    pknn: float,
) -> float:
    return nt / (14000.0 * dil) * vsdil * (per(ph_f, pknn) - per(ph_s, pknn))


def m_hpo4(
    ph_f: float,
    ph_s: float,
    pt: float,
    dil: float,
    vsdil: float,
    pkpp: float,
) -> float:
    return pt / (31000.0 * dil) * vsdil * (per(ph_f, pkpp) - per(ph_s, pkpp))
