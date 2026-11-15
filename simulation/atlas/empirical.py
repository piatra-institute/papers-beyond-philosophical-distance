"""Higher-order structure in two public datasets of positions on philosophical questions."""
from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from .maxent import (LOG2, MaxEnt, benjamini_hochberg, connected_information_3, g2_batch,
                     g2_three_way, independent_fit, ipf_pairwise, marginal3_from_full, states,
                     table3, three_way_log_interaction)

POSITION_ITEMS = ["driving", "light", "phenomenal", "higher-order", "explanatorygap", "biology",
                  "evolution", "machines", "zombie", "mary", "transition"]
ENTITY_ITEMS = ["entities_amoeba", "entities_baby", "entities_bat", "entities_dog", "entities_fish",
                "entities_monkey", "entities_octopus", "entities_other_people", "entities_rock",
                "entities_tree", "entities_thermostat", "entities_worm", "entities_yourself"]
LAMBDAS = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
CELL_LABELS = ["".join(map(str, c)) for c in itertools.product((0, 1), repeat=3)]


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
    out |= {"folds": folds,
            "gain_pairwise_over_independent_bits": float(d21.mean()),
            "gain_third_order_over_pairwise_bits": float(d32.mean()),
            "gain_third_order_over_pairwise_se_bits": float(d32.std(ddof=1) / math.sqrt(folds)),
            "folds_third_order_better": int((d32 > 0).sum())}
    return out


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
        "strata": strata,
        "bootstrap": {"replicates": boot, "seed": seed},
    }


# --- survey of consciousness researchers ----------------------------------------------

def load_items(path: Path, cols: list[str], cut: int = 2) -> tuple[np.ndarray, list[str]]:
    """Binary endorsement (answer 1 or 2 on the five-point scale), complete cases."""
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if all(r[c] != "" for c in cols):
                rows.append([1.0 if int(r[c]) <= cut else 0.0 for c in cols])
    return np.array(rows), cols


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
        out[f"gain_order{a}_over_order{b}_se_bits"] = float(d.std(ddof=1) / math.sqrt(len(d)))
        out[f"folds_order{a}_better"] = int((d > 0).sum())
    for o in orders:
        lams = [r[f"lambda_{o}"] for r in recs]
        out[f"lambda_mode_order{o}"] = max(set(lams), key=lams.count)
    return out


def triplet_tests(X: np.ndarray, names: list[str], base: MaxEnt, nsim: int, seed: int, q: float = 0.05) -> dict:
    """Third-order test of each triplet against the fitted pairwise model of all items.

    The null distribution of G^2 for a triplet is simulated from the triplet's marginal
    under the fitted pairwise model (which can carry small third-order terms induced by
    marginalising over the other items).
    """
    rng = np.random.default_rng(seed)
    n = X.shape[1]
    N = len(X)
    p_full = base.probabilities()
    rows = []
    for i, j, k in itertools.combinations(range(n), 3):
        t = table3(X, i, j, k)
        g2_obs, p_asym = g2_three_way(t)
        marg = marginal3_from_full(p_full, n, i, j, k)
        sims = rng.multinomial(N, marg.ravel(), size=nsim).reshape(nsim, 2, 2, 2)
        g2_null = g2_batch(sims)
        p_boot = (1 + int((g2_null >= g2_obs - 1e-12).sum())) / (nsim + 1)
        ci = connected_information_3(t)
        rows.append({"items": [names[i], names[j], names[k]], "g2": g2_obs, "p_chi2": p_asym, "p_boot": p_boot,
                     "third_order_bits": ci["third_order_bits"], "third_order_fraction": ci["third_order_fraction"],
                     "null_mean_third_order_bits": float(g2_null.mean() / (2 * N * LOG2)),
                     "log_interaction": three_way_log_interaction(t, pseudo=0.5)})
    pb = np.array([r["p_boot"] for r in rows])
    disc = benjamini_hochberg(pb, q)
    top = sorted(rows, key=lambda r: r["p_boot"])[:3]
    excess = np.array([r["third_order_bits"] - r["null_mean_third_order_bits"] for r in rows])
    return {"triplets": len(rows), "simulations_per_triplet": nsim, "fdr_q": q,
            "discoveries": int(disc.sum()), "min_p_boot": float(pb.min()),
            "p_boot_below_005": int((pb < 0.05).sum()),
            "median_third_order_bits": float(np.median([r["third_order_bits"] for r in rows])),
            "median_null_third_order_bits": float(np.median([r["null_mean_third_order_bits"] for r in rows])),
            "median_excess_third_order_bits": float(np.median(excess)),
            "max_third_order_bits": float(max(r["third_order_bits"] for r in rows)),
            "median_third_order_fraction": float(np.median([r["third_order_fraction"] for r in rows])),
            "top": top}


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
    """Detection power of the triplet test for planted third-order structure.

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
            g2 = g2_batch(sims)
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
                     nsim_avoid: int, power_ref_bits: float | None, extreme: float = 0.05) -> dict:
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
    with open(path, newline="", encoding="utf-8") as fh:
        total_rows = sum(1 for _ in csv.DictReader(fh))
    out = {"label": label, "respondents_total": total_rows, "respondents_complete": N, "items": names,
           "items_count": len(names), "dropped_extreme_items": dropped, "dropped_count": len(dropped),
           "extreme_threshold": extreme,
           "endorsement_rates": {c: float(r) for c, r in zip(names, rates)},
           "cross_validation": cv,
           "base_pairwise_lambda": lam,
           "triplet_tests": triplet_tests(X, names, base, nsim_triplets, seed + 2),
           "joint_avoidances": avoidance_test(X, base, nsim_avoid, seed + 3)}
    if power_ref_bits is not None:
        n_tests = math.comb(n, 3)
        out["power"] = power_analysis(base, N, n_triplets=8, steps=24, nsim=4000, seed=seed + 4,
                                      n_tests=n_tests, reference_bits=power_ref_bits)
    return out


def design_power(sample_sizes: dict[str, int], triplets: dict[str, int], reference_bits: float,
                 simulated: dict | None = None) -> dict:
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
    return out
