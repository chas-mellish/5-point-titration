from __future__ import annotations

import argparse
import json
import sys

from .constants import THESIS_DEFAULTS
from .core import run_titration
from .models import TitrationInput


def _build_parser() -> argparse.ArgumentParser:
    d = THESIS_DEFAULTS
    p = argparse.ArgumentParser(
        prog="titration",
        description="5-point titration analysis for anaerobic digester water chemistry",
    )
    p.add_argument("--ph0", type=float, default=d.ph0, help=f"initial sample pH (default: {d.ph0})")
    p.add_argument("--ph1", type=float, default=d.ph1, help=f"pH at titration point 1 (default: {d.ph1})")
    p.add_argument("--ph2", type=float, default=d.ph2, help=f"pH at titration point 2 (default: {d.ph2})")
    p.add_argument("--ph3", type=float, default=d.ph3, help=f"pH at titration point 3 (default: {d.ph3})")
    p.add_argument("--ph4", type=float, default=d.ph4, help=f"pH at titration point 4 (default: {d.ph4})")
    p.add_argument("--vx1", type=float, default=d.vx1, help=f"titrant volume at point 1, mL (default: {d.vx1})")
    p.add_argument("--vx2", type=float, default=d.vx2, help=f"titrant volume at point 2, mL (default: {d.vx2})")
    p.add_argument("--vx3", type=float, default=d.vx3, help=f"titrant volume at point 3, mL (default: {d.vx3})")
    p.add_argument("--vx4", type=float, default=d.vx4, help=f"titrant volume at point 4, mL (default: {d.vx4})")
    p.add_argument("--titrant-normality", type=float, default=d.titrant_normality, help=f"titrant normality, eq/L (default: {d.titrant_normality})")
    p.add_argument("--sample-volume-undiluted", type=float, default=d.sample_volume_undiluted, help=f"undiluted sample volume, mL (default: {d.sample_volume_undiluted})")
    p.add_argument("--sample-volume-diluted", type=float, default=d.sample_volume_diluted, help=f"diluted sample volume, mL (default: {d.sample_volume_diluted})")
    p.add_argument("--temperature", type=float, default=d.temperature, help=f"sample temperature, deg C (default: {d.temperature})")
    p.add_argument("--tds", type=float, default=d.tds, help=f"total dissolved solids, mg/L (default: {d.tds})")
    p.add_argument("--inorganic-nitrogen", type=float, default=d.inorganic_nitrogen, help=f"inorganic nitrogen, mg/L as N (default: {d.inorganic_nitrogen})")
    p.add_argument("--inorganic-phosphorus", type=float, default=d.inorganic_phosphorus, help=f"inorganic phosphorus, mg/L as P (default: {d.inorganic_phosphorus})")
    p.add_argument("--json", action="store_true", help="Output results as JSON instead of formatted text")
    return p


def _format_results(result) -> str:
    lines = [
        "",
        "===== 5-Point Titration Results =====",
        "",
        f"  H2CO3* alkalinity:          {result.h2co3_alkalinity:10.2f}  mg/L as CaCO3",
        f"  Short-chain fatty acids:    {max(result.scfa_concentration, 0):10.2f}  mg/L as acetic acid",
        f"  Systematic pH error:        {result.systematic_ph_error:10.4f}",
        "",
    ]

    if result.convergence_status == "converged":
        lines.append("  Convergence: OK")
    elif result.convergence_status == "exceeded_max_iterations":
        lines.append("  WARNING: Maximum iterations exceeded — result may be inaccurate")
    elif result.convergence_status == "ratio_too_high":
        lines.append("  WARNING: SCFA/alkalinity ratio too high — iteration skipped")

    lines.append("")
    return "\n".join(lines)


def _format_json(result, inp: TitrationInput) -> str:
    """Convert results to JSON-friendly dictionary"""
    return json.dumps({
        "input": {
            "ph_values": [inp.ph0, inp.ph1, inp.ph2, inp.ph3, inp.ph4],
            "volumes_mL": [inp.vx1, inp.vx2, inp.vx3, inp.vx4],
            "titrant_normality_eq_L": inp.titrant_normality,
            "sample_volume_undiluted_mL": inp.sample_volume_undiluted,
            "sample_volume_diluted_mL": inp.sample_volume_diluted,
            "temperature_C": inp.temperature,
            "tds_mg_L": inp.tds,
            "inorganic_nitrogen_mg_L_as_N": inp.inorganic_nitrogen,
            "inorganic_phosphorus_mg_L_as_P": inp.inorganic_phosphorus,
        },
        "results": {
            "h2co3_alkalinity_mg_L_as_CaCO3": round(result.h2co3_alkalinity, 2),
            "scfa_concentration_mg_L_as_acetic": round(max(result.scfa_concentration, 0), 2),
            "systematic_ph_error": round(result.systematic_ph_error, 4),
            "ct_values": [round(result.ct_values[0], 2), round(result.ct_values[1], 2)],
            "at_values": [round(result.at_values[0], 2), round(result.at_values[1], 2)],
            "convergence_status": result.convergence_status,
        },
        "metadata": {
            "tool": "5-point-titration",
            "version": "0.1.0",
            "pascal_source": "Moosbrugger 1991 doctoral thesis",
        },
    }, indent=2)


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    inp = TitrationInput(
        ph0=args.ph0,
        ph1=args.ph1,
        ph2=args.ph2,
        ph3=args.ph3,
        ph4=args.ph4,
        vx1=args.vx1,
        vx2=args.vx2,
        vx3=args.vx3,
        vx4=args.vx4,
        titrant_normality=args.titrant_normality,
        sample_volume_undiluted=args.sample_volume_undiluted,
        sample_volume_diluted=args.sample_volume_diluted,
        temperature=args.temperature,
        tds=args.tds,
        inorganic_nitrogen=args.inorganic_nitrogen,
        inorganic_phosphorus=args.inorganic_phosphorus,
    )

    result = run_titration(inp)
    
    if args.json:
        sys.stdout.write(_format_json(result, inp))
    else:
        sys.stdout.write(_format_results(result))


if __name__ == "__main__":
    main()
