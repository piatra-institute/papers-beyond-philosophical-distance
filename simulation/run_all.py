"""Recorded run: tests, exhaustive verification of obstruction orders, data analysis, figures.

Writes output/results.json and ../paper/figures/{orders,data}.pdf. Reads only the derived
files in data/. Run fetch_data.py first: it rebuilds the per-respondent survey file, which is
not committed, and checks every derived file against its pinned SHA-256.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
import numpy as np
import scipy

from atlas import empirical as E
from atlas import figures as Fg
from atlas import lpbounds, theory

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "output"
FIG = HERE.parent / "paper" / "figures"
SEED = 20261115


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def to_jsonable(x):
    if isinstance(x, dict):
        return {str(k): to_jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [to_jsonable(v) for v in x]
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, float) and (x != x):
        return None
    return x


def main() -> int:
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    tests = subprocess.run([sys.executable, "-m", "pytest", "-q", "tests"], cwd=HERE,
                           capture_output=True, text=True)
    tests_passed = tests.returncode == 0
    manifest = json.loads((DATA / "MANIFEST.json").read_text())
    inputs = {name: sha256(DATA / name) for name in manifest["derived"]}
    inputs_ok = all(inputs[n] == manifest["derived"][n]["sha256"] for n in inputs)

    res: dict = {"seed": SEED}
    res["theory"] = theory.run(planar_m=range(4, 17), planar_sets=100)
    res["theory"]["lp_exhaustive"] = lpbounds.exhaustive_comparison()
    res["trolley"] = E.trolley_analysis(DATA / "trolley_first_sessions.csv", boot=2000, seed=SEED)
    ref_bits = res["trolley"]["information"]["third_order_bits"]
    ref_lam = res["trolley"]["log_interaction"]
    discovery = {"targets": [ref_bits, 0.02, 0.04, 0.06, 0.08, 0.10], "n_triplets": 4, "reps": 25,
                 "null_reps": 100}
    res["consciousness"] = {
        "positions": E.battery_analysis(DATA / "consciousness_items.csv", E.POSITION_ITEMS, "positions",
                                        seed=SEED, nsim_triplets=9999, nsim_avoid=500, power_ref_bits=ref_bits,
                                        reference_log_interaction=ref_lam, discovery=discovery),
        "entities": E.battery_analysis(DATA / "consciousness_items.csv", E.ENTITY_ITEMS, "entities",
                                       seed=SEED + 100, nsim_triplets=9999, nsim_avoid=500, power_ref_bits=None,
                                       reference_log_interaction=ref_lam),
    }
    pos = res["consciousness"]["positions"]
    res["design"] = E.design_power(
        {"survey_of_consciousness_researchers": pos["respondents_complete"], "philpapers_2020_full_sample": 1785},
        {"survey_of_consciousness_researchers": pos["power"]["tests_for_bonferroni"],
         "philpapers_2020_full_sample": 161700},
        reference_bits=ref_bits, simulated=pos["power"], discovery=pos["discovery_power"])
    Fg.figure_orders(res, FIG / "orders.pdf")
    Fg.figure_data(res, FIG / "data.pdf")
    th = res["theory"]
    pos_lp = res["consciousness"]["positions"]["lp_bounds"]
    ent_lp = res["consciousness"]["entities"]["lp_bounds"]
    checks = {
        "tests_passed": tests_passed,
        "inputs_match_manifest": inputs_ok,
        "boolean_formulas_match": th["boolean"]["all_formulas_match"],
        "boolean_methods_agree": th["boolean"]["all_methods_agree"],
        "bundling_invariant": th["bundling"]["invariant"],
        "witnesses_verified": th["witnesses"]["all_verified"],
        "single_issue_and_products_match": th["single_issue_and_products"]["all_match"],
        "planar_order_equals_empty_polygon_all": th["planar"]["order_equals_largest_empty_polygon"] == th["planar"]["sets_total"],
        "convex_position_order_equals_m": th["planar"]["convex_position_order_equals_m"],
        "background_logic_all_agree": th["background_logic"]["all_agree"],
        "pairwise_hull_all_agree": th["pairwise_hull"]["all_agree"],
        "parity_pairs_equal": th["parity_pairs"]["all_equal"],
        "lp_exhaustive_all_agree": th["lp_exhaustive"]["all_agree"],
        "lp_data_certificates_verified": all(
            b["certificates_verified"] == b["certificates_attempted"] and b["facial_set_agrees_with_support_lp"]
            and b["hull_is_majority_closure"] for b in (pos_lp, ent_lp, res["trolley"]["lp_bounds"])),
        "monte_carlo_resolution_sufficient": all(
            res["consciousness"][k]["triplet_tests"]["smallest_attainable_p"]
            < res["consciousness"][k]["triplet_tests"]["first_bh_threshold"] for k in ("positions", "entities")),
        "planted_margins_preserved": res["consciousness"]["positions"]["discovery_power"]["max_margin_error"] < 1e-6,
    }
    res["checks"] = checks
    res["environment"] = {"python": platform.python_version(), "numpy": np.__version__,
                          "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
                          "platform": platform.platform()}
    res["inputs"] = inputs
    res["tests_output_tail"] = tests.stdout.strip().splitlines()[-3:]
    res["runtime_seconds"] = round(time.time() - t0, 1)
    res["execution"] = {"status": "ok" if all(checks.values()) else "check_failed"}
    (OUT / "results.json").write_text(json.dumps(to_jsonable(res), indent=1) + "\n", encoding="utf-8")
    print(json.dumps(checks, indent=1))
    print("runtime", res["runtime_seconds"], "s")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
