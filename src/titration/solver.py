from __future__ import annotations

from .chemistry import (
    calculate_pk_constants,
    d_h2co3_alk,
    d_hac_alk,
    m_h2o,
    m_hpo4,
    m_nh3,
    per_co3,
    per_hco3,
)
from .constants import ACETIC_ACID_FACTOR, CACO3_FACTOR
from .models import TitrationInput, TitrationResult


def _to_mg_per_l(molar_qty, factor, vsdil, dil):
    return molar_qty / vsdil * factor * dil


def run_solver(inp: TitrationInput) -> TitrationResult:
    if inp.sample_volume_undiluted <= 0:
        raise ValueError("sample_volume_undiluted must be positive")
    if inp.sample_volume_diluted <= 0:
        raise ValueError("sample_volume_diluted must be positive")

    dil = inp.sample_volume_diluted / inp.sample_volume_undiluted
    vsdil = inp.sample_volume_diluted
    ca = inp.titrant_normality
    nt = inp.inorganic_nitrogen
    pt = inp.inorganic_phosphorus

    pks = calculate_pk_constants(inp.temperature, inp.tds, dil)
    pk11 = pks["pk11"]
    pk22 = pks["pk22"]
    pkaa = pks["pkaa"]
    pknn = pks["pknn"]
    pkpp = pks["pkpp"]
    logf1 = pks["logf1"]

    # Working copies of pH values -- mutated during iteration
    ph0 = inp.ph0
    ph1 = inp.ph1
    ph2 = inp.ph2
    ph3 = inp.ph3
    ph4 = inp.ph4

    if ph3 == ph4:
        raise ValueError("ph3 and ph4 must differ (division by d_hac_alk(ph3, ph4) is zero)")
    if ph1 == ph4:
        raise ValueError("ph1 and ph4 must differ (division by d_hac_alk(ph1, ph4) is zero)")

    vx1 = inp.vx1
    vx2 = inp.vx2
    vx3 = inp.vx3
    vx4 = inp.vx4

    def _atctcalculation():
        def _compute_ct_pair(ph_a, ph_b, vx_a, vx_b):
            ratio = d_hac_alk(ph_a, ph_b, pkaa) / d_hac_alk(ph3, ph4, pkaa)
            a = (
                (vx_b - vx_a) * ca
                - m_h2o(vx_a, vx_b, ph_a, ph_b, vsdil, logf1)
                - m_nh3(ph_a, ph_b, nt, dil, vsdil, pknn)
                - m_hpo4(ph_a, ph_b, pt, dil, vsdil, pkpp)
                + ratio
                * (
                    m_h2o(vx3, vx4, ph3, ph4, vsdil, logf1)
                    + m_nh3(ph3, ph4, nt, dil, vsdil, pknn)
                    + m_hpo4(ph3, ph4, pt, dil, vsdil, pkpp)
                    - (vx4 - vx3) * ca
                )
            )
            b = d_h2co3_alk(ph_a, ph_b, pk11, pk22) - ratio * d_h2co3_alk(
                ph3, ph4, pk11, pk22
            )
            return a, b

        a1, b1 = _compute_ct_pair(ph1, ph2, vx1, vx2)
        a2, b2 = _compute_ct_pair(ph1, ph4, vx1, vx4)

        if b1 == 0.0:
            raise ValueError("b1 is zero — degenerate pH combination in Ct1 computation")
        if b2 == 0.0:
            raise ValueError("b2 is zero — degenerate pH combination in Ct2 computation")

        m_ct1 = a1 / b1
        ct1 = _to_mg_per_l(m_ct1, CACO3_FACTOR, vsdil, dil)
        ct2 = _to_mg_per_l(a2 / b2, CACO3_FACTOR, vsdil, dil)
        ct_comp = ct1 - ct2
        return m_ct1, ct1, ct2, ct_comp

    # Initial calculation
    del_ph = 0.0
    ph_corr = 0.0
    counter = 0

    m_ct1, ct1, ct2, ct_comp = _atctcalculation()

    def _compute_m_at(ph_a, ph_b, vx_a, vx_b, m_ct1_val):
        return (1.0 / d_hac_alk(ph_a, ph_b, pkaa)) * (
            (vx_b - vx_a) * ca
            - m_ct1_val * d_h2co3_alk(ph_a, ph_b, pk11, pk22)
            - m_nh3(ph_a, ph_b, nt, dil, vsdil, pknn)
            - m_hpo4(ph_a, ph_b, pt, dil, vsdil, pkpp)
            - m_h2o(vx_a, vx_b, ph_a, ph_b, vsdil, logf1)
        )

    # Initial At1 and ratio for direction check
    m_at1 = _compute_m_at(ph3, ph4, vx3, vx4, m_ct1)
    at1_initial = _to_mg_per_l(m_at1, ACETIC_ACID_FACTOR, vsdil, dil)
    ct_at_ratio = at1_initial / ct1

    # Determine correction direction
    if ct_at_ratio > 0.5:
        convergence_status = "ratio_too_high"
    elif ct_comp == 0.0:
        convergence_status = "converged"
    else:
        ph_corr = -0.01 if ct_comp > 0.0 else 0.01
        initial_sign = 1.0 if ct_comp > 0.0 else -1.0
        while not (ct_comp * initial_sign < 0.0 or counter > 19):
            del_ph += ph_corr
            ph0 += ph_corr
            ph1 += ph_corr
            ph2 += ph_corr
            ph3 += ph_corr
            ph4 += ph_corr
            counter += 1
            m_ct1, ct1, ct2, ct_comp = _atctcalculation()
        convergence_status = (
            "exceeded_max_iterations" if counter > 19 else "converged"
        )

    # Final At calculations using converged MCt1 and adjusted pH values
    m_at1 = _compute_m_at(ph3, ph4, vx3, vx4, m_ct1)
    at1 = _to_mg_per_l(m_at1, ACETIC_ACID_FACTOR, vsdil, dil)

    m_at2 = _compute_m_at(ph1, ph4, vx1, vx4, m_ct1)
    at2 = _to_mg_per_l(m_at2, ACETIC_ACID_FACTOR, vsdil, dil)

    # H2CO3* alkalinity of the undiluted sample
    m_h2co3_alk_sam = m_ct1 * (per_hco3(ph0, pk11, pk22) + 2.0 * per_co3(ph0, pk11, pk22))
    h2co3_alk_sam = (
        _to_mg_per_l(m_h2co3_alk_sam, CACO3_FACTOR, vsdil, dil)
        + (10.0 ** (ph0 - 14.0) - 10.0 ** (-ph0) / 10.0 ** logf1) * CACO3_FACTOR
    )

    # Clamp SCFA to zero if negative, matching Pascal output procedure
    scfa = max(at1, 0.0)

    return TitrationResult(
        h2co3_alkalinity=h2co3_alk_sam,
        scfa_concentration=scfa,
        systematic_ph_error=del_ph,
        convergence_status=convergence_status,
        ct_values=(ct1, ct2),
        at_values=(at1, at2),
    )
