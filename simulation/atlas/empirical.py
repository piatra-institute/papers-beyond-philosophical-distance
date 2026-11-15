"""Higher-order structure in two public datasets of positions on philosophical questions."""
from __future__ import annotations

import csv
import itertools
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.stats import chi2, norm
from scipy.stats import t as student_t

from . import lpbounds
from .maxent import (LOG2, MaxEnt, benjamini_hochberg, connected_information_3, g2_fast,
                     g2_three_way, independent_fit, ipf_pairwise, marginal3_from_full, table3,
                     three_way_log_interaction)

POSITION_ITEMS = ["driving", "light", "phenomenal", "higher-order", "explanatorygap", "biology",
                  "evolution", "machines", "zombie", "mary", "transition"]
ENTITY_ITEMS = ["entities_amoeba", "entities_baby", "entities_bat", "entities_dog", "entities_fish",
                "entities_monkey", "entities_octopus", "entities_other_people", "entities_rock",
                "entities_tree", "entities_thermostat", "entities_worm", "entities_yourself"]
LAMBDAS = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
CELL_LABELS = ["".join(map(str, c)) for c in itertools.product((0, 1), repeat=3)]


# --- uncertainty of cross-validated differences ---------------------------------------

def corrected_cv_summary(d: np.ndarray, test_sizes: np.ndarray, n_total: int) -> dict:
    """Mean of per-fold differences with the corrected resampled variance.

    Fold scores share respondents across training sets and across repetitions, so the
    naive variance s^2 / J understates the uncertainty. Nadeau and Bengio's correction
    replaces 1/J by 1/J + n2/n1, with n2 and n1 the test and training sizes; the interval
    uses Student's t with J - 1 degrees of freedom. The fold standard deviation is
    reported as a descriptive quantity.
    """
    d = np.asarray(d, dtype=float)
    J = len(d)
    n2 = float(np.mean(test_sizes))
    n1 = float(n_total - n2)
    var = float(d.var(ddof=1))
    se = math.sqrt((1.0 / J + n2 / n1) * var)
    tq = float(student_t.ppf(0.975, J - 1))
    return {"mean": float(d.mean()), "se_corrected": se, "ci95_corrected": [float(d.mean() - tq * se), float(d.mean() + tq * se)],
            "sd_folds": math.sqrt(var), "folds": J, "test_over_train": n2 / n1, "folds_positive": int((d > 0).sum())}


# --- trolley dilemmas ----------------------------------------------------------------

def load_trolley(path: Path, field: str = "country_code") -> dict[str, np.ndarray]:
    """2x2x2 tables (Switch, Loop, Footbridge; 1 = sacrifice accepted) by a stratifying field."""
    tables: dict[str, np.ndarray] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            t = tables.setdefault(row[field], np.zeros((2, 2, 2)))
            t[int(row["switch"]), int(row["loop"]), int(row["footbridge"])] += int(row["count"])
    return tables


def stratified_third_order(path: Path, field: str, min_n: int, lam_pooled: float,
                           valid=lambda key: True) -> dict:
    """Third-order structure within strata of a design or sample variable."""
    tables = {k: t for k, t in load_trolley(path, field).items() if valid(k) and t.sum() >= min_n}
    rows = []
    for key, t in sorted(tables.items()):
        ci = connected_information_3(t)
        g2, p = g2_three_way(t)
        rows.append({"stratum": key, "n": int(t.sum()), "third_order_bits": ci["third_order_bits"],
                     "third_order_fraction": ci["third_order_fraction"],
                     "log_interaction": three_way_log_interaction(t), "g2": g2, "p": p})
    k = len(rows)
    n_tot = sum(r["n"] for r in rows)
    sub = sum(tables.values()) if tables else np.zeros((2, 2, 2))
    return {"field": field, "min_participants": min_n, "strata": k, "participants": n_tot,
            "same_sign_as_pooled": int(sum(1 for r in rows if np.sign(r["log_interaction"]) == np.sign(lam_pooled))),
            "p_below_005": int(sum(1 for r in rows if r["p"] < 0.05)),
            "p_below_bonferroni": int(sum(1 for r in rows if r["p"] < 0.05 / max(k, 1))),
            "median_third_order_fraction": float(np.median([r["third_order_fraction"] for r in rows])) if rows else float("nan"),
            "within_stratum_third_order_bits_weighted": float(sum(r["n"] * r["third_order_bits"] for r in rows) / n_tot) if n_tot else float("nan"),
            "combined_third_order_bits": connected_information_3(sub)["third_order_bits"] if n_tot else float("nan"),
            "rows": rows}


def _is_country(code: str) -> bool:
    return len(code) == 2 and code.isalpha() and code.isupper()


def _cv_three_items(table: np.ndarray, folds: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    idx = np.repeat(np.arange(8), table.ravel().astype(int))
    rng.shuffle(idx)
    parts = np.array_split(idx, folds)
    scores = {"independent": [], "pairwise": [], "saturated": []}
    sizes = []
    for f in range(folds):
        test = parts[f]
        train = np.concatenate([parts[g] for g in range(folds) if g != f])
        tr = np.bincount(train, minlength=8).reshape(2, 2, 2).astype(float)
        models = {"independent": independent_fit(tr), "pairwise": ipf_pairwise(tr),
                  "saturated": (tr + 0.5) / (tr.sum() + 4.0)}
        for k, q in models.items():
            scores[k].append(float(np.mean(np.log2(q.ravel()[test]))))
        sizes.append(len(test))
    out = {k: float(np.average(v, weights=sizes)) for k, v in scores.items()}
    d32 = np.array(scores["saturated"]) - np.array(scores["pairwise"])
    d21 = np.array(scores["pairwise"]) - np.array(scores["independent"])
    total = int(table.sum())
    out |= {"folds": folds,
            "gain_pairwise_over_independent_bits": float(d21.mean()),
            "gain_third_order_over_pairwise_bits": float(d32.mean()),
            "gain_third_order_over_pairwise": corrected_cv_summary(d32, np.array(sizes), total),
            "gain_pairwise_over_independent": corrected_cv_summary(d21, np.array(sizes), total),
            "folds_third_order_better": int((d32 > 0).sum())}
    return out


def _pattern_counts_trolley(table: np.ndarray) -> np.ndarray:
    """Counts indexed by pattern integers with bit 0 = Switch, bit 1 = Loop, bit 2 = Footbridge."""
    counts = np.zeros(8, dtype=np.int64)
    for a, b, c in itertools.product((0, 1), repeat=3):
        counts[a | (b << 1) | (c << 2)] = int(table[a, b, c])
    return counts


def trolley_analysis(path: Path, boot: int = 2000, seed: int = 20261115, min_country: int = 1000) -> dict:
    tables = load_trolley(path)
    pooled = sum(tables.values())
    N = int(pooled.sum())
    ci = connected_information_3(pooled)
    g2, p = g2_three_way(pooled)
    fit = ipf_pairwise(pooled)
    obs = pooled / N
    ratios = obs / fit
    rng = np.random.default_rng(seed)
    samples = rng.multinomial(N, obs.ravel(), size=boot).reshape(boot, 2, 2, 2).astype(float)
    fits = ipf_pairwise(samples)
    probs = samples / N
    with np.errstate(divide="ignore", invalid="ignore"):
        h3 = -np.where(probs > 0, probs * np.log2(np.where(probs > 0, probs, 1.0)), 0.0).sum(axis=(1, 2, 3))
        h2 = -np.where(fits > 0, fits * np.log2(np.where(fits > 0, fits, 1.0)), 0.0).sum(axis=(1, 2, 3))
    m0 = probs.sum(axis=(2, 3)); m1 = probs.sum(axis=(1, 3)); m2 = probs.sum(axis=(1, 2))
    ind = m0[:, :, None, None] * m1[:, None, :, None] * m2[:, None, None, :]
    h1 = -(ind * np.log2(ind)).sum(axis=(1, 2, 3))
    i3_b = h2 - h3
    frac_b = i3_b / (h1 - h3)
    ratio_b = (probs / fits).reshape(boot, 8)
    lam_b = np.array([three_way_log_interaction(s) for s in samples])
    cells = []
    for c, lab in enumerate(CELL_LABELS):
        a, b, d = (int(ch) for ch in lab)
        cells.append({"pattern_switch_loop_footbridge": lab, "count": int(pooled[a, b, d]),
                      "observed": float(obs[a, b, d]), "pairwise_prediction": float(fit[a, b, d]),
                      "ratio": float(ratios[a, b, d]),
                      "ratio_ci95": [float(np.percentile(ratio_b[:, c], 2.5)), float(np.percentile(ratio_b[:, c], 97.5))]})
    lam_pooled = three_way_log_interaction(pooled)
    strata = {
        "country": stratified_third_order(path, "country_code", min_country, lam_pooled, valid=_is_country),
        "presentation_order": stratified_third_order(path, "order", min_country, lam_pooled),
        "sacrifice_side": stratified_third_order(path, "sacrifice_side_SLF", min_country, lam_pooled),
        "device": stratified_third_order(path, "template", min_country, lam_pooled),
        "after_moral_machine": stratified_third_order(path, "after_moral_machine", min_country, lam_pooled),
    }
    return {
        "participants": N,
        "country_codes_valid": int(sum(1 for c in tables if _is_country(c))),
        "participants_without_valid_country": int(sum(t.sum() for c, t in tables.items() if not _is_country(c))),
        "acceptance": {"switch": float(obs[1].sum()), "loop": float(obs[:, 1].sum()),
                       "footbridge": float(obs[:, :, 1].sum())},
        "information": ci,
        "third_order_bits_ci95": [float(np.percentile(i3_b, 2.5)), float(np.percentile(i3_b, 97.5))],
        "third_order_fraction_ci95": [float(np.percentile(frac_b, 2.5)), float(np.percentile(frac_b, 97.5))],
        "g2": g2, "g2_p": p,
        "log_interaction": lam_pooled,
        "log_interaction_ci95": [float(np.percentile(lam_b, 2.5)), float(np.percentile(lam_b, 97.5))],
        "cells": cells,
        "cross_validation": _cv_three_items(pooled, folds=10, seed=seed),
        "lp_bounds": {k: v for k, v in lpbounds.analyse_counts(_pattern_counts_trolley(pooled), 3).items()
                      if k != "excluded_patterns"},
        "strata": strata,
        "bootstrap": {"replicates": boot, "seed": seed},
    }


# --- survey of consciousness researchers ----------------------------------------------

def eligible(row: dict) -> bool:
    """The screening of Francken et al. (2022): self-identified consciousness researchers
    (consciousness = 1) who had not filled in the survey before (before != 1)."""
    return row["consciousness"] == "1" and row["before"] != "1"


def screening_counts(path: Path) -> dict:
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    first = [r for r in rows if r["before"] != "1"]
    return {"submissions": len(rows),
            "repeat_submissions": sum(1 for r in rows if r["before"] == "1"),
            "not_researcher": sum(1 for r in first if r["consciousness"] == "2"),
            "researcher_status_missing": sum(1 for r in first if r["consciousness"] == ""),
            "eligible": sum(1 for r in rows if eligible(r)),
            "eligible_attended_assc": sum(1 for r in rows if eligible(r) and r["assc"] == "1")}


def load_items(path: Path, cols: list[str], cut: int = 2) -> tuple[np.ndarray, list[str]]:
    """Binary endorsement (answer 1 or 2 on the five-point scale), eligible complete cases."""
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if eligible(r) and all(r[c] != "" for c in cols):
                rows.append([1.0 if int(r[c]) <= cut else 0.0 for c in cols])
    return np.array(rows), cols


def five_point_answers(path: Path, col: str, cols: list[str]) -> int:
    """Eligible complete respondents who answered 2 or 4 on an item whose online forms had
    three answers coded 1, 3 and 5 (possible only on the five-point form)."""
    with open(path, newline="", encoding="utf-8") as fh:
        return sum(1 for r in csv.DictReader(fh)
                   if eligible(r) and all(r[c] != "" for c in cols) and r[col] in ("2", "4"))


def select_lambda(X: np.ndarray, order: int, lambdas: list[float], folds: int, rng: np.random.Generator) -> float:
    n = X.shape[1]
    perm = rng.permutation(len(X))
    parts = np.array_split(perm, folds)
    best, best_lam = -np.inf, lambdas[0]
    for lam in lambdas:
        tot = 0.0
        for f in range(folds):
            te = parts[f]
            tr = np.concatenate([parts[g] for g in range(folds) if g != f])
            tot += MaxEnt(n, order).fit(X[tr], lam).mean_loglik_bits(X[te]) * len(te)
        if tot > best:
            best, best_lam = tot, lam
    return best_lam


def cv_compare(X: np.ndarray, seed: int, reps: int = 5, outer: int = 10, inner: int = 5,
               lambdas: list[float] = LAMBDAS, orders=(1, 2, 3)) -> dict:
    """Nested cross-validation: penalty chosen on the training part of each outer fold."""
    rng = np.random.default_rng(seed)
    n = X.shape[1]
    recs = []
    for _ in range(reps):
        perm = rng.permutation(len(X))
        parts = np.array_split(perm, outer)
        for f in range(outer):
            te = parts[f]
            tr = np.concatenate([parts[g] for g in range(outer) if g != f])
            rec = {"n_test": len(te)}
            for o in orders:
                lam = select_lambda(X[tr], o, lambdas, inner, rng)
                rec[o] = MaxEnt(n, o).fit(X[tr], lam).mean_loglik_bits(X[te])
                rec[f"lambda_{o}"] = lam
            recs.append(rec)
    w = np.array([r["n_test"] for r in recs], dtype=float)
    out = {"reps": reps, "outer_folds": outer, "inner_folds": inner, "lambdas": lambdas,
           "folds_total": len(recs)}
    for o in orders:
        out[f"heldout_bits_order{o}"] = float(np.average([r[o] for r in recs], weights=w))
    for a, b in ((2, 1), (3, 2)):
        d = np.array([r[a] - r[b] for r in recs])
        out[f"gain_order{a}_over_order{b}_bits"] = float(d.mean())
        out[f"gain_order{a}_over_order{b}"] = corrected_cv_summary(d, w, len(X))
        out[f"folds_order{a}_better"] = int((d > 0).sum())
    for o in orders:
        lams = [r[f"lambda_{o}"] for r in recs]
        out[f"lambda_mode_order{o}"] = max(set(lams), key=lams.count)
    return out


def triplet_pvalues(X: np.ndarray, base: MaxEnt, nsim: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Monte Carlo p-values of the third-order test for every triplet of items.

    For each triplet, G^2 of the observed 2x2x2 table against its no-three-way fit is
    compared with nsim tables drawn from the triplet's marginal under the fitted pairwise
    model of all items (which can carry small third-order terms induced by marginalising
    over the other items). The smallest attainable p-value is 1 / (nsim + 1).
    """
    n = X.shape[1]
    N = len(X)
    p_full = base.probabilities()
    trips = list(itertools.combinations(range(n), 3))
    g2_obs = g2_fast(np.stack([table3(X, i, j, k) for i, j, k in trips]))
    pv = np.empty(len(trips))
    null_mean = np.empty(len(trips))
    for t_idx, (i, j, k) in enumerate(trips):
        marg = marginal3_from_full(p_full, n, i, j, k)
        sims = rng.multinomial(N, marg.ravel() / marg.sum(), size=nsim).reshape(nsim, 2, 2, 2)
        g2_null = g2_fast(sims)
        pv[t_idx] = (1 + int((g2_null >= g2_obs[t_idx] - 1e-9).sum())) / (nsim + 1)
        null_mean[t_idx] = g2_null.mean()
    return pv, g2_obs, null_mean


def triplet_tests(X: np.ndarray, names: list[str], base: MaxEnt, nsim: int, seed: int, q: float = 0.05,
                  reference_log_interaction: float | None = None) -> dict:
    """Third-order test of each triplet, with Benjamini-Hochberg control and simultaneous
    intervals for the three-way log interaction."""
    rng = np.random.default_rng(seed)
    n = X.shape[1]
    N = len(X)
    pb, g2_obs, g2_null_mean = triplet_pvalues(X, base, nsim, rng)
    trips = list(itertools.combinations(range(n), 3))
    m = len(trips)
    zsim = float(norm.isf(0.025 / m))
    rows = []
    for t_idx, (i, j, k) in enumerate(trips):
        t = table3(X, i, j, k)
        ci = connected_information_3(t)
        lam = three_way_log_interaction(t, pseudo=0.5)
        se = math.sqrt(float(np.sum(1.0 / (t + 0.5))))
        null_mean_bits = float(g2_null_mean[t_idx] / (2 * N * LOG2))
        rows.append({"items": [names[i], names[j], names[k]], "g2": float(g2_obs[t_idx]),
                     "p_chi2": float(chi2.sf(g2_obs[t_idx], 1)), "p_boot": float(pb[t_idx]),
                     "third_order_bits": ci["third_order_bits"], "third_order_fraction": ci["third_order_fraction"],
                     "null_mean_third_order_bits": null_mean_bits,
                     "log_interaction": lam, "log_interaction_se": se,
                     "log_interaction_simultaneous_ci95": [lam - zsim * se, lam + zsim * se]})
    disc = benjamini_hochberg(pb, q)
    top = sorted(rows, key=lambda r: r["p_boot"])[:3]
    excess = np.array([r["third_order_bits"] - r["null_mean_third_order_bits"] for r in rows])
    half = np.array([zsim * r["log_interaction_se"] for r in rows])
    out = {"triplets": m, "simulations_per_triplet": nsim, "fdr_q": q,
           "smallest_attainable_p": 1.0 / (nsim + 1), "first_bh_threshold": q / m,
           "discoveries": int(disc.sum()), "min_p_boot": float(pb.min()),
           "p_boot_below_005": int((pb < 0.05).sum()),
           "expected_below_005_under_null": 0.05 * m,
           "median_third_order_bits": float(np.median([r["third_order_bits"] for r in rows])),
           "median_null_third_order_bits": float(np.median([r["null_mean_third_order_bits"] for r in rows])),
           "median_excess_third_order_bits": float(np.median(excess)),
           "max_third_order_bits": float(max(r["third_order_bits"] for r in rows)),
           "median_third_order_fraction": float(np.median([r["third_order_fraction"] for r in rows])),
           "simultaneous_z": zsim,
           "simultaneous_halfwidth_median": float(np.median(half)),
           "simultaneous_halfwidth_min": float(half.min()),
           "max_abs_log_interaction": float(max(abs(r["log_interaction"]) for r in rows)),
           "top": top}
    if reference_log_interaction is not None:
        ref = abs(reference_log_interaction)
        out["reference_log_interaction"] = ref
        out["triplets_excluding_reference_size"] = int(sum(
            1 for r in rows if -ref < r["log_interaction_simultaneous_ci95"][0] and r["log_interaction_simultaneous_ci95"][1] < ref))
        out["triplets_interval_contains_zero"] = int(sum(
            1 for r in rows if r["log_interaction_simultaneous_ci95"][0] <= 0 <= r["log_interaction_simultaneous_ci95"][1]))
    return out


# --- discovery power of the full screening procedure ---------------------------------------

def planted_distribution(base: MaxEnt, trip: tuple[int, int, int], beta: float) -> tuple[np.ndarray, float]:
    """Distribution over all answer patterns with the one- and two-way margins of `base`
    and an added three-way term beta * x_i x_j x_k (the I-projection of the tilted model
    onto those margins, fitted by convex optimisation)."""
    i, j, k = trip
    S, F = base.S, base.F
    mu = base.probabilities() @ F
    off = beta * S[:, i] * S[:, j] * S[:, k]

    def obj(th: np.ndarray) -> tuple[float, np.ndarray]:
        e = F @ th + off
        mx = e.max()
        w = np.exp(e - mx)
        z = w.sum()
        return float(np.log(z) + mx - th @ mu), (w / z) @ F - mu

    res = minimize(obj, base.theta.copy(), jac=True, method="L-BFGS-B",
                   options={"maxiter": 20000, "gtol": 1e-13, "ftol": 1e-16})
    e = F @ res.x + off
    p = np.exp(e - e.max())
    p /= p.sum()
    return p, float(np.max(np.abs(p @ F - mu)))


def _planted_bits(base: MaxEnt, trip, beta) -> float:
    p, _ = planted_distribution(base, trip, beta)
    return connected_information_3(marginal3_from_full(p, base.n, *trip))["third_order_bits"]


def _screen_once(task: tuple) -> tuple[int, int, float]:
    """Simulate one dataset from a planted model and rerun the screening procedure:
    penalty chosen by inner cross-validation, pairwise fit, Monte Carlo p-values for all
    triplets, Benjamini-Hochberg at q. Returns (planted triplet discovered, other
    discoveries, p-value of the planted triplet)."""
    p, n, N, trip_index, nsim, q, seed = task
    rng = np.random.default_rng(seed)
    S = np.array(list(itertools.product((0, 1), repeat=n)), dtype=np.float64)
    X = S[rng.choice(len(p), size=N, p=p)]
    lam = select_lambda(X, 2, LAMBDAS, 5, rng)
    base = MaxEnt(n, 2).fit(X, lam)
    pv, _, _ = triplet_pvalues(X, base, nsim, rng)
    disc = benjamini_hochberg(pv, q)
    hit = int(disc[trip_index]) if trip_index >= 0 else 0
    return hit, int(disc.sum()) - hit, float(pv[trip_index]) if trip_index >= 0 else float("nan")


def discovery_power(base: MaxEnt, N: int, targets: list[float], n_triplets: int, reps: int, nsim: int,
                    seed: int, q: float = 0.05, null_reps: int = 100, workers: int | None = None) -> dict:
    """Probability that the full screening procedure discovers a planted third-order term.

    For randomly chosen triplets, a three-way term is added to the fitted pairwise model
    of the battery while every one- and two-way margin is held fixed, with strength set
    so that the triplet's own third-order information equals each target. Each simulated
    dataset of N respondents is analysed exactly as the observed data are. The null level
    (no planted term) estimates the probability of any false discovery.
    """
    rng = np.random.default_rng(seed)
    n = base.n
    trips = list(itertools.combinations(range(n), 3))
    chosen = [int(c) for c in rng.choice(len(trips), size=n_triplets, replace=False)]
    p0 = base.probabilities()
    betas = np.linspace(0.0, 8.0, 161)
    plan = []
    for ti in chosen:
        trip = trips[ti]
        bits = np.array([_planted_bits(base, trip, b) for b in betas])
        start = int(np.argmin(bits))
        row = {"triplet": [int(v) for v in trip], "levels": []}
        for target in targets:
            idx = next((a for a in range(start, len(betas) - 1) if bits[a] < target <= bits[a + 1]), None)
            if idx is None:
                row["levels"].append(None)
                continue
            beta = brentq(lambda b: _planted_bits(base, trip, b) - target, betas[idx], betas[idx + 1], xtol=1e-10)
            pp, err = planted_distribution(base, trip, beta)
            row["levels"].append({"target_bits": target, "beta": float(beta), "margin_error": err, "p": pp})
        plan.append((ti, row))
    seeds = np.random.SeedSequence(seed).generate_state(len(targets) * n_triplets * reps + null_reps)
    tasks, keys = [], []
    s_i = 0
    for r in range(null_reps):
        tasks.append((p0, n, N, -1, nsim, q, int(seeds[s_i]))); keys.append(("null", None)); s_i += 1
    for ti, row in plan:
        for li, lev in enumerate(row["levels"]):
            for r in range(reps):
                if lev is not None:
                    tasks.append((lev["p"], n, N, ti, nsim, q, int(seeds[s_i]))); keys.append((li, ti))
                s_i += 1
    workers = workers or max(1, min(8, (os.cpu_count() or 2) - 1))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        results = list(ex.map(_screen_once, tasks, chunksize=4))
    null_any = [res[1] > 0 for key, res in zip(keys, results) if key[0] == "null"]
    levels = []
    for li, target in enumerate(targets):
        rs = [res for key, res in zip(keys, results) if key[0] == li]
        hits = np.array([r_[0] for r_ in rs])
        k_, n_ = int(hits.sum()), len(hits)
        zc = 1.959964
        centre = (k_ + zc * zc / 2) / (n_ + zc * zc)
        half = zc * math.sqrt(k_ * (n_ - k_) / n_ + zc * zc / 4) / (n_ + zc * zc) if n_ else float("nan")
        levels.append({"target_bits": target, "datasets": n_, "discovered": k_,
                       "power": k_ / n_ if n_ else float("nan"),
                       "power_ci95": [max(0.0, centre - half), min(1.0, centre + half)],
                       "mean_other_discoveries": float(np.mean([r_[1] for r_ in rs])) if rs else float("nan"),
                       "median_planted_p": float(np.median([r_[2] for r_ in rs])) if rs else float("nan")})
    pw = [lv["power"] for lv in levels]
    bits80 = None
    for a in range(len(targets) - 1):
        if pw[a] < 0.8 <= pw[a + 1]:
            bits80 = targets[a] + (0.8 - pw[a]) * (targets[a + 1] - targets[a]) / (pw[a + 1] - pw[a])
            break
    return {"sample_size": N, "triplets_planted": n_triplets, "replicates_per_triplet_and_level": reps,
            "simulations_per_triplet_test": nsim, "fdr_q": q, "workers_independent_seeds": True,
            "datasets_total": len(tasks),
            "null_datasets": len(null_any), "null_any_discovery": float(np.mean(null_any)) if null_any else None,
            "levels": levels, "bits_for_80pct_discovery": bits80,
            "max_margin_error": max(lev["margin_error"] for _, row in plan for lev in row["levels"] if lev),
            "planted": [{"triplet": row["triplet"],
                         "betas": [lev["beta"] if lev else None for lev in row["levels"]]} for _, row in plan]}


def joint_avoidances(X: np.ndarray) -> int:
    """Triplet answer patterns never observed although each of their pairs is observed."""
    n = X.shape[1]
    total = 0
    for i, j, k in itertools.combinations(range(n), 3):
        t = table3(X, i, j, k)
        pij, pik, pjk = t.sum(axis=2), t.sum(axis=1), t.sum(axis=0)
        for a, b, c in itertools.product((0, 1), repeat=3):
            if t[a, b, c] == 0 and pij[a, b] > 0 and pik[a, c] > 0 and pjk[b, c] > 0:
                total += 1
    return total


def avoidance_test(X: np.ndarray, base: MaxEnt, nsim: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    obs = joint_avoidances(X)
    sims = np.array([joint_avoidances(base.sample(len(X), rng)) for _ in range(nsim)])
    return {"observed": obs, "simulations": nsim, "null_mean": float(sims.mean()),
            "null_ci95": [float(np.percentile(sims, 2.5)), float(np.percentile(sims, 97.5))],
            "p_upper": float((1 + (sims >= obs).sum()) / (nsim + 1))}


def power_analysis(base: MaxEnt, N: int, n_triplets: int, steps: int, nsim: int, seed: int,
                   n_tests: int, reference_bits: float) -> dict:
    """Power of the test for one pre-specified triplet, without multiplicity correction.

    This is the power of a likelihood-ratio test of a single triplet chosen in advance,
    against the asymptotic chi-square cut-off at alpha = 0.05 (and, for reference, at the
    Bonferroni level over all triplets); the power of the screening procedure actually
    applied to the survey is estimated by discovery_power.

    For each randomly chosen triplet, q* is the pairwise (no three-way interaction) fit
    of the triplet's marginal under the fitted pairwise model of all items. The planted
    table is q* + eps * s, where s(a, b, c) = (-1)^(a + b + c): this direction leaves
    every one- and two-way margin unchanged, so the planted third-order information is
    exactly KL(q || q*) and grows monotonically with eps up to the positivity limit.
    """
    rng = np.random.default_rng(seed)
    n = base.n
    trip_all = list(itertools.combinations(range(n), 3))
    chosen = [trip_all[i] for i in rng.choice(len(trip_all), size=n_triplets, replace=False)]
    p_full = base.probabilities()
    sgn = np.array([[[(-1.0) ** (a + b + c) for c in (0, 1)] for b in (0, 1)] for a in (0, 1)])
    crit05 = float(chi2.isf(0.05, 1))
    critb = float(chi2.isf(0.05 / n_tests, 1))
    rows = []
    for (i, j, k) in chosen:
        qstar = ipf_pairwise(marginal3_from_full(p_full, n, i, j, k))
        eps_max = 0.98 * float(qstar[sgn < 0].min())
        curve = []
        for step in range(steps + 1):
            eps = eps_max * step / steps
            q = qstar + eps * sgn
            i3 = float(np.sum(q * np.log2(q / qstar)))
            sims = rng.multinomial(N, q.ravel(), size=nsim).reshape(nsim, 2, 2, 2)
            g2 = g2_fast(sims)
            curve.append({"eps": eps, "third_order_bits": i3,
                          "power_alpha_005": float((g2 > crit05).mean()),
                          "power_alpha_bonferroni": float((g2 > critb).mean())})
        rows.append({"triplet": [int(i), int(j), int(k)], "eps_max": eps_max, "curve": curve})

    def bits_at(curve, key, target):
        xs = [c["third_order_bits"] for c in curve]
        ys = [c[key] for c in curve]
        for a in range(len(ys) - 1):
            if ys[a] < target <= ys[a + 1]:
                t = (target - ys[a]) / (ys[a + 1] - ys[a])
                return xs[a] + t * (xs[a + 1] - xs[a])
        return float("nan")

    def power_at(curve, key, x):
        xs = [c["third_order_bits"] for c in curve]
        ys = [c[key] for c in curve]
        return float(np.interp(x, xs, ys))

    b80 = [bits_at(r["curve"], "power_alpha_005", 0.8) for r in rows]
    b80b = [bits_at(r["curve"], "power_alpha_bonferroni", 0.8) for r in rows]
    pref = [power_at(r["curve"], "power_alpha_005", reference_bits) for r in rows]
    size = [r["curve"][0]["power_alpha_005"] for r in rows]
    finite = lambda v: [x for x in v if x == x]
    return {"triplets": n_triplets, "steps": steps, "simulations": nsim, "sample_size": N,
            "tests_for_bonferroni": n_tests,
            "bits_for_80pct_power_alpha005_median": float(np.median(finite(b80))) if finite(b80) else None,
            "bits_for_80pct_power_alpha005_range": [min(finite(b80)), max(finite(b80))] if finite(b80) else None,
            "bits_for_80pct_power_bonferroni_median": float(np.median(finite(b80b))) if finite(b80b) else None,
            "bits_for_80pct_power_bonferroni_range": [min(finite(b80b)), max(finite(b80b))] if finite(b80b) else None,
            "triplets_reaching_80pct_bonferroni": len(finite(b80b)),
            "reference_bits": reference_bits,
            "power_at_reference_median": float(np.median(pref)),
            "power_at_reference_range": [float(min(pref)), float(max(pref))],
            "size_at_zero_mean": float(np.mean(size)),
            "rows": rows}


def battery_analysis(path: Path, cols: list[str], label: str, seed: int, nsim_triplets: int,
                     nsim_avoid: int, power_ref_bits: float | None, extreme: float = 0.05,
                     reference_log_interaction: float | None = None, discovery: dict | None = None) -> dict:
    X, names = load_items(path, cols)
    rates = X.mean(axis=0)
    keep = [c for c, r in zip(names, rates) if extreme <= r <= 1 - extreme]
    dropped = [c for c in names if c not in keep]
    if dropped:
        X, names = load_items(path, keep)
        rates = X.mean(axis=0)
    N, n = X.shape
    cv = cv_compare(X, seed=seed)
    rng = np.random.default_rng(seed + 1)
    lam = select_lambda(X, 2, LAMBDAS, 5, rng)
    base = MaxEnt(n, 2).fit(X, lam)
    codes = (X.astype(np.int64) * (1 << np.arange(n))).sum(axis=1)
    lp = lpbounds.analyse_counts(np.bincount(codes, minlength=1 << n), n)
    out = {"label": label, "screening": screening_counts(path), "respondents_complete": N, "items": names,
           "items_count": len(names), "dropped_extreme_items": dropped, "dropped_count": len(dropped),
           "extreme_threshold": extreme,
           "endorsement_rates": {c: float(r) for c, r in zip(names, rates)},
           "cross_validation": cv,
           "base_pairwise_lambda": lam,
           "lp_bounds": {k: v for k, v in lp.items() if k != "excluded_patterns"},
           "triplet_tests": triplet_tests(X, names, base, nsim_triplets, seed + 2,
                                          reference_log_interaction=reference_log_interaction),
           "joint_avoidances": avoidance_test(X, base, nsim_avoid, seed + 3)}
    if "transition" in names:
        out["transition_answers_2_or_4"] = five_point_answers(path, "transition", names)
    if power_ref_bits is not None:
        n_tests = math.comb(n, 3)
        out["power"] = power_analysis(base, N, n_triplets=8, steps=24, nsim=4000, seed=seed + 4,
                                      n_tests=n_tests, reference_bits=power_ref_bits)
    if discovery is not None:
        out["discovery_power"] = discovery_power(base, N, seed=seed + 5, nsim=nsim_triplets, **discovery)
    return out


def design_power(sample_sizes: dict[str, int], triplets: dict[str, int], reference_bits: float,
                 simulated: dict | None = None, discovery: dict | None = None) -> dict:
    """Analytic power of the triplet test (noncentral chi-square, 1 df).

    Under a third-order effect of I bits, G^2 is approximately noncentral chi-square with
    one degree of freedom and noncentrality 2 N ln(2) I. Returns, for each design, the
    effect detectable with 80% power at alpha = 0.05 (one pre-specified triplet) and with
    a Bonferroni correction over all triplets, and the power at the reference effect.
    """
    from scipy.optimize import brentq
    from scipy.stats import ncx2

    def lam80(crit: float) -> float:
        return brentq(lambda lam: ncx2.sf(crit, 1, lam) - 0.8, 1e-6, 1e4)

    out = {}
    for name, N in sample_sizes.items():
        crit05 = float(chi2.isf(0.05, 1))
        critb = float(chi2.isf(0.05 / triplets[name], 1))
        scale = 2.0 * N * LOG2
        lam_ref = scale * reference_bits
        out[name] = {"sample_size": N, "triplets": triplets[name],
                     "bits_for_80pct_power_alpha005": lam80(crit05) / scale,
                     "bits_for_80pct_power_bonferroni": lam80(critb) / scale,
                     "power_at_reference_alpha005": float(ncx2.sf(crit05, 1, lam_ref)),
                     "power_at_reference_bonferroni": float(ncx2.sf(critb, 1, lam_ref))}
    out["reference_bits"] = reference_bits
    if simulated is not None:
        out["check_against_simulation"] = {
            "analytic_alpha005": out["survey_of_consciousness_researchers"]["bits_for_80pct_power_alpha005"],
            "simulated_alpha005_median": simulated["bits_for_80pct_power_alpha005_median"],
            "analytic_bonferroni": out["survey_of_consciousness_researchers"]["bits_for_80pct_power_bonferroni"],
            "simulated_bonferroni_median": simulated["bits_for_80pct_power_bonferroni_median"]}
        if discovery is not None:
            out["check_against_simulation"]["simulated_screening_80pct_discovery"] = discovery["bits_for_80pct_discovery"]
    return out
