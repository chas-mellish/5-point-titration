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


def run_solver(inp: TitrationInput) -> TitrationResult:
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

    vx1 = inp.vx1
    vx2 = inp.vx2
    vx3 = inp.vx3
    vx4 = inp.vx4

    def _atctcalculation():
        ratio_12_34 = d_hac_alk(ph1, ph2, pkaa) / d_hac_alk(ph3, ph4, pkaa)
        a1 = (
            (vx2 - vx1) * ca
            - m_h2o(vx1, vx2, ph1, ph2, vsdil, logf1)
            - m_nh3(ph1, ph2, nt, dil, vsdil, pknn)
            - m_hpo4(ph1, ph2, pt, dil, vsdil, pkpp)
            + ratio_12_34
            * (
                m_h2o(vx3, vx4, ph3, ph4, vsdil, logf1)
                + m_nh3(ph3, ph4, nt, dil, vsdil, pknn)
                + m_hpo4(ph3, ph4, pt, dil, vsdil, pkpp)
                - (vx4 - vx3) * ca
            )
        )
        b1 = d_h2co3_alk(ph1, ph2, pk11, pk22) - ratio_12_34 * d_h2co3_alk(
            ph3, ph4, pk11, pk22
        )

        ratio_14_34 = d_hac_alk(ph1, ph4, pkaa) / d_hac_alk(ph3, ph4, pkaa)
        a2 = (
            (vx4 - vx1) * ca
            - m_h2o(vx1, vx4, ph1, ph4, vsdil, logf1)
            - m_nh3(ph1, ph4, nt, dil, vsdil, pknn)
            - m_hpo4(ph1, ph4, pt, dil, vsdil, pkpp)
            + ratio_14_34
            * (
                m_h2o(vx3, vx4, ph3, ph4, vsdil, logf1)
                + m_nh3(ph3, ph4, nt, dil, vsdil, pknn)
                + m_hpo4(ph3, ph4, pt, dil, vsdil, pkpp)
                - (vx4 - vx3) * ca
            )
        )
        b2 = d_h2co3_alk(ph1, ph4, pk11, pk22) - ratio_14_34 * d_h2co3_alk(
            ph3, ph4, pk11, pk22
        )

        m_ct1 = a1 / b1
        ct1 = m_ct1 / vsdil * CACO3_FACTOR * dil
        ct2 = (a2 / b2) / vsdil * CACO3_FACTOR * dil
        ct_comp = ct1 - ct2
        return m_ct1, ct1, ct2, ct_comp

    # Initial calculation
    del_ph = 0.0
    ph_corr = 0.0
    counter = 0

    m_ct1, ct1, ct2, ct_comp = _atctcalculation()

    # Initial At1 and ratio for direction check
    m_at1 = (1.0 / d_hac_alk(ph3, ph4, pkaa)) * (
        (vx4 - vx3) * ca
        - m_ct1 * d_h2co3_alk(ph3, ph4, pk11, pk22)
        - m_nh3(ph3, ph4, nt, dil, vsdil, pknn)
        - m_hpo4(ph3, ph4, pt, dil, vsdil, pkpp)
        - m_h2o(vx3, vx4, ph3, ph4, vsdil, logf1)
    )
    at1_initial = m_at1 / vsdil * ACETIC_ACID_FACTOR * dil
    ct_at_ratio = at1_initial / ct1

    # Determine correction direction
    if ct_at_ratio > 0.5:
        convergence_status = "ratio_too_high"
    elif ct_comp == 0.0:
        convergence_status = "converged"
    elif ct_comp > 0.0:
        # Case 'a': need negative pH correction
        ph_corr = -0.01
        while not (ct_comp < 0.0 or counter > 19):
            del_ph -= 0.01
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
    else:
        # Case 'b': need positive pH correction
        ph_corr = 0.01
        while not (ct_comp > 0.0 or counter > 19):
            del_ph += 0.01
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
    m_at1 = (1.0 / d_hac_alk(ph3, ph4, pkaa)) * (
        (vx4 - vx3) * ca
        - m_ct1 * d_h2co3_alk(ph3, ph4, pk11, pk22)
        - m_nh3(ph3, ph4, nt, dil, vsdil, pknn)
        - m_hpo4(ph3, ph4, pt, dil, vsdil, pkpp)
        - m_h2o(vx3, vx4, ph3, ph4, vsdil, logf1)
    )
    at1 = m_at1 / vsdil * ACETIC_ACID_FACTOR * dil

    m_at2 = (1.0 / d_hac_alk(ph1, ph4, pkaa)) * (
        (vx4 - vx1) * ca
        - m_ct1 * d_h2co3_alk(ph1, ph4, pk11, pk22)
        - m_nh3(ph1, ph4, nt, dil, vsdil, pknn)
        - m_hpo4(ph1, ph4, pt, dil, vsdil, pkpp)
        - m_h2o(vx1, vx4, ph1, ph4, vsdil, logf1)
    )
    at2 = m_at2 / vsdil * ACETIC_ACID_FACTOR * dil

    # H2CO3* alkalinity of the undiluted sample
    m_h2co3_alk_sam = m_ct1 * (per_hco3(ph0, pk11, pk22) + 2.0 * per_co3(ph0, pk11, pk22))
    h2co3_alk_sam = (
        m_h2co3_alk_sam / vsdil * dil * CACO3_FACTOR
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
