"""Recorded run: tests, exhaustive verification of obstruction orders, data analysis, figures.

Writes output/results.json and ../paper/figures/{orders,data}.pdf. Reads only the frozen
files in data/ (see fetch_data.py, which is not run here).
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
from atlas import theory

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
    res["trolley"] = E.trolley_analysis(DATA / "trolley_first_sessions.csv", boot=2000, seed=SEED)
    ref_bits = res["trolley"]["information"]["third_order_bits"]
    res["consciousness"] = {
        "positions": E.battery_analysis(DATA / "consciousness_items.csv", E.POSITION_ITEMS, "positions",
                                        seed=SEED, nsim_triplets=2000, nsim_avoid=500, power_ref_bits=ref_bits),
        "entities": E.battery_analysis(DATA / "consciousness_items.csv", E.ENTITY_ITEMS, "entities",
                                       seed=SEED + 100, nsim_triplets=2000, nsim_avoid=500, power_ref_bits=None),
    }
    pos = res["consciousness"]["positions"]
    res["design"] = E.design_power(
        {"survey_of_consciousness_researchers": pos["respondents_complete"], "philpapers_2020_full_sample": 1785},
        {"survey_of_consciousness_researchers": pos["power"]["tests_for_bonferroni"],
         "philpapers_2020_full_sample": 161700},
        reference_bits=ref_bits, simulated=pos["power"])
    Fg.figure_orders(res, FIG / "orders.pdf")
    Fg.figure_data(res, FIG / "data.pdf")
    th = res["theory"]
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
