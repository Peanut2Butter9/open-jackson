"""Recompute every number and regenerate the slides' LaTeX snippets.

    python -m openjackson                     # full run, about two minutes
    python -m openjackson --caps 10 20 --reps 2 --horizon 20000   # quick check
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import ctmc, hospital, paper, tex
from .simulate import simulate

STATIONS = ["Registrasi", "Dokter", "Laboratorium"]


def thousands(x) -> str:
    """12345 -> 12\\,345 (a thin space, so it reads the same with . or , as decimal mark)."""
    return f"{int(x):,}".replace(",", r"\,")


def sci(x) -> str:
    """0.00025 -> 2.5\\times10^{-4}, for use in math mode."""
    mantissa, exponent = f"{x:.1e}".split("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def hospital_results(caps, reps, horizon) -> dict:
    net = hospital.network()
    exact = net.mean_number()
    convergence = []
    for cap in caps:
        pi = ctmc.stationary_distribution(net, cap)
        small = np.ix_(*[range(6)] * 3)  # states with at most 5 patients per station
        product = np.array([[[net.probability((a, b, c)) for c in range(6)] for b in range(6)] for a in range(6)])
        convergence.append(
            {
                "cap": cap,
                "states": int(pi.size),
                "mean_number": ctmc.mean_number(pi).tolist(),
                "max_rel_error_small_states": float(np.abs(pi[small] / product - 1).max()),
            }
        )
    runs = [simulate(net, horizon, seed=s) for s in range(reps)]
    sim_number = np.array([r.mean_number for r in runs])
    sim_time = np.array([r.mean_time_in_system for r in runs])
    return {
        "symbolic": {
            "arrival_rates": [str(x) for x in hospital.arrival_rates()],
            "utilisation": [str(x) for x in hospital.utilisation()],
            "cases": [
                {
                    "case": str(case),
                    "product_form_satisfies_derived_equation": hospital.product_form_holds(case),
                    "original_equation_matches": hospital.same_equation(
                        hospital.ORIGINAL_EQUATIONS[case], hospital.balance_equation(case)
                    ),
                    "product_form_satisfies_original_equation": hospital.product_form_holds(
                        case, hospital.ORIGINAL_EQUATIONS[case]
                    ),
                }
                for case in hospital.CASES
            ],
        },
        "example": {
            "rates": {str(k): v for k, v in hospital.EXAMPLE.items()},
            "arrival_rates": net.arrival_rates().tolist(),
            "utilisation": net.utilisation().tolist(),
            "mean_number": exact.tolist(),
            "time_per_visit": net.time_per_visit().tolist(),
            "mean_number_in_system": net.mean_number_in_system(),
            "mean_time_in_system": net.mean_time_in_system(),
            "sum_of_time_per_visit": float(net.time_per_visit().sum()),
            "ctmc": convergence,
            "simulation": {
                "replications": reps,
                "horizon_minutes": horizon,
                "warmup_fraction": 0.1,
                "mean_number": sim_number.mean(axis=0).tolist(),
                "mean_number_se": (sim_number.std(axis=0, ddof=1) / np.sqrt(reps)).tolist(),
                "mean_time_in_system": float(sim_time.mean()),
                "mean_time_in_system_se": float(sim_time.std(ddof=1) / np.sqrt(reps)),
                "patients": int(sum(r.patients for r in runs)),
            },
        },
    }


def write_tex(out: Path, paper_res: dict, hosp: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    f = tex.number
    (out / "transitions.tex").write_text(tex.transition_rows())
    for k, case in enumerate(hospital.CASES, start=1):
        (out / f"diagram-{k}.tex").write_text(tex.diagram(case))
        (out / f"balance-{k}.tex").write_text(tex.balance(case))
        (out / f"case-{k}.tex").write_text(f"${tex.condition(case)}$")

    lam, mu, rho = paper_res["lambda"], paper_res["mu"], paper_res["rho"]
    rows = [
        f"{i} & {a} & {d} & {f(lam[i - 1])} & {f(mu[i - 1])} & {f(paper.REPORTED_RHO[i - 1])} & {f(rho[i - 1])} \\\\"
        for i, a, d in zip(range(1, 7), paper.ARRIVALS, paper.DEPARTURES)
    ]
    (out / "paper-rates.tex").write_text("\n".join(rows) + "\n")

    ls, norm = paper_res["routing_least_squares"], paper_res["routing_normalised"]
    rows = [
        f"$P_{{{name[1:]}}}$ & {f(ls[name])} & {f(norm[name], 7)} & {f(rep, 7)} \\\\"
        for name, rep in paper.REPORTED_ROUTING.items()
    ]
    (out / "paper-routing.tex").write_text("\n".join(rows) + "\n")

    per = paper_res["mean_number_per_station"]
    t_lam, t_rho = paper_res["traffic_arrival_rates"], paper_res["traffic_utilisation"]
    bold = lambda x: rf"\textbf{{{f(x)}}}" if x >= 1 else f(x)  # noqa: E731
    rows = [
        f"{i} & {f(lam[i - 1])} & {f(per[i - 1], 2)} & {round(per[i - 1])} & {f(per[i - 1] / lam[i - 1], 2)}"
        f" & {f(t_lam[i - 1])} & {bold(t_rho[i - 1])} \\\\"
        for i in range(1, 7)
    ]
    (out / "paper-stations.tex").write_text("\n".join(rows) + "\n")

    ex = hosp["example"]
    sim, last = ex["simulation"], ex["ctmc"][-1]
    rows = [
        f"{i} {name} & {f(ex['arrival_rates'][i - 1])} & {f(ex['utilisation'][i - 1])} & {f(ex['mean_number'][i - 1], 3)}"
        f" & {f(last['mean_number'][i - 1], 3)} & {f(sim['mean_number'][i - 1], 3)} $\\pm$ {f(sim['mean_number_se'][i - 1], 3)} \\\\"
        for i, name in enumerate(STATIONS, start=1)
    ]
    rows.append(
        f"\\midrule Total $L$ & & & {f(ex['mean_number_in_system'], 3)} & {f(sum(last['mean_number']), 3)}"
        f" & {f(sum(sim['mean_number']), 3)} \\\\"
    )
    (out / "example.tex").write_text("\n".join(rows) + "\n")

    rows = [
        f"{c['cap']} & {thousands(c['states'])} & {' & '.join(f(x, 4) for x in c['mean_number'])}"
        f" & {f(sum(c['mean_number']), 4)} & ${sci(c['max_rel_error_small_states'])}$ \\\\"
        for c in ex["ctmc"]
    ]
    (out / "ctmc.tex").write_text("\n".join(rows) + "\n")

    rho_sym = hospital.utilisation()
    macros = {
        "PaperL": f(paper_res["mean_number"], 2),
        "PaperWLittle": f(paper_res["mean_time_little"], 2),
        "PaperRhoFive": f(t_rho[4], 3),
        "PaperPinvDiff": sci(paper_res["pseudo_inverse_max_abs_diff"]),
        "RhoOne": sp.latex(rho_sym[0]),
        "RhoTwo": sp.latex(rho_sym[1]),
        "RhoThree": sp.latex(rho_sym[2]),
        "ExLamOne": f(float(hospital.EXAMPLE[hospital.lam1])),
        "ExMuOne": f(float(hospital.EXAMPLE[hospital.mu1])),
        "ExMuTwo": f(float(hospital.EXAMPLE[hospital.mu2])),
        "ExMuThree": f(float(hospital.EXAMPLE[hospital.mu3])),
        "ExL": f(ex["mean_number_in_system"], 3),
        "ExW": f(ex["mean_time_in_system"], 2),
        "ExWSum": f(ex["sum_of_time_per_visit"], 2),
        "ExWSim": f(sim["mean_time_in_system"], 2),
        "ExWSimSe": f(sim["mean_time_in_system_se"], 2),
        "ExSimReps": sim["replications"],
        "ExSimHorizon": thousands(sim["horizon_minutes"]),
        "ExSimPatients": thousands(sim["patients"]),
        "CtmcCap": last["cap"],
    }
    (out / "numbers.tex").write_text(tex.macros(macros))


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--tex", type=Path, default=Path("slides/generated"))
    parser.add_argument("--caps", type=int, nargs="+", default=[10, 20, 30, 40, 50])
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--horizon", type=int, default=200_000, help="simulated minutes per replication")
    args = parser.parse_args(argv)

    paper_res = paper.reproduce()
    hosp = hospital_results(args.caps, args.reps, args.horizon)
    args.results.mkdir(parents=True, exist_ok=True)
    reported = {
        "lambda": paper.REPORTED_LAMBDA.tolist(),
        "mu": paper.REPORTED_MU.tolist(),
        "rho": paper.REPORTED_RHO.tolist(),
        "routing": paper.REPORTED_ROUTING,
        "mean_number": paper.REPORTED_MEAN_NUMBER,
        "mean_time": paper.REPORTED_MEAN_TIME,
    }
    (args.results / "paper.json").write_text(json.dumps({"reported": reported, "reproduced": paper_res}, indent=2) + "\n")
    (args.results / "hospital.json").write_text(json.dumps(hosp, indent=2) + "\n")
    write_tex(args.tex, paper_res, hosp)
    print(f"wrote {args.results}/paper.json, {args.results}/hospital.json and {args.tex}/*.tex")


if __name__ == "__main__":
    main()
