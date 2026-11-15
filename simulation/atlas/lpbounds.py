"""What numerical pairwise marginals exclude: linear-programming bounds on answer patterns.

Patterns are integers p in [0, 2^n) whose bit i is the answer to item i (the convention
of atlas.theory). For binary items the pairwise marginals of a distribution p are its
moments E[1], E[x_i] and E[x_i x_j], collected as Phi^T p, where Phi is the 2^n x m
matrix of monomials of degree at most two. For a pattern z,

    max p(z)  subject to  p >= 0 and Phi^T p = mu                         (primal)
  = min mu . y  subject to  Phi y >= e_z                                   (dual)

A dual-feasible y is a pairwise function f = Phi y with f(z) >= 1 and f >= 0 elsewhere,
so every distribution with marginals mu satisfies p(z) <= E_p[f] = mu . y. A zero
optimum is certified by an exact rational y with mu . y = 0, checked here in integer
arithmetic. The patterns with a positive optimum form the facial set of mu, the support
of every distribution in the relative interior of the set of distributions with these
marginals (Fienberg and Rinaldo 2012). It contains the observed support S, depends on
the data only through S, and lies inside the pairwise hull H(S), the majority closure
of S. The facial set is also computed from S directly by a single linear program, as an
independent check of the per-pattern bounds.
"""
from __future__ import annotations

import itertools
import math
import random
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

from .theory import majority_closure, pairwise_hull

ZERO_TOL = 1e-9


def moment_matrix(n: int) -> np.ndarray:
    """Rows: patterns; columns: 1, x_i (i < n), x_i x_j (i < j). Integer entries."""
    P = np.arange(1 << n)
    bits = (P[:, None] >> np.arange(n)[None, :]) & 1
    cols = [np.ones(1 << n, dtype=np.int64)] + [bits[:, i] for i in range(n)]
    cols += [bits[:, i] * bits[:, j] for i, j in itertools.combinations(range(n), 2)]
    return np.stack(cols, axis=1).astype(np.int64)


def max_probability(Phi: np.ndarray, mu: np.ndarray, z: int) -> float:
    """Largest probability of pattern z among distributions with moments mu (sum mu_0 = 1)."""
    c = np.zeros(Phi.shape[0])
    c[z] = -1.0
    res = linprog(c, A_eq=Phi.T.astype(float), b_eq=np.asarray(mu, dtype=float), bounds=(0, None),
                  method="highs")
    if res.status != 0:
        raise RuntimeError(f"primal LP failed for pattern {z}: {res.message}")
    return float(-res.fun)


def _exact_check(Phi: np.ndarray, mu_counts: np.ndarray, z: int, y: list[Fraction]) -> tuple[bool, int]:
    L = math.lcm(*[q.denominator for q in y])
    yi = [int(q * L) for q in y]
    big = max(abs(v) for v in yi) if yi else 0
    if big < 2 ** 40:
        f = Phi @ np.array(yi, dtype=np.int64)
        f = [int(v) for v in f]
    else:
        f = [sum(int(a) * b for a, b in zip(row, yi)) for row in Phi]
    feasible = all(v >= (L if x == z else 0) for x, v in enumerate(f))
    value = sum(int(a) * b for a, b in zip(mu_counts, yi))
    return feasible and value == 0, L


def certificate(Phi: np.ndarray, mu_counts: np.ndarray, z: int, max_den: int = 10 ** 6) -> dict:
    """Dual optimum for pattern z and an exact check of the certificate it yields.

    mu_counts are the integer moments of the observed counts. The floating dual optimum
    is rounded to nearby rationals; the certificate is accepted only if the rounded
    pairwise function is nonnegative on every pattern, at least one on z, and has
    expectation exactly zero under the observed marginals.
    """
    e = np.zeros(Phi.shape[0])
    e[z] = 1.0
    scale = float(mu_counts[0])
    res = linprog(np.asarray(mu_counts, dtype=float) / scale, A_ub=-Phi.astype(float), b_ub=-e,
                  bounds=(None, None), method="highs")
    if res.status != 0:
        return {"pattern": int(z), "dual_value": None, "verified": False, "denominator": None}
    for snap in (1e-9, 1e-7):
        y = [Fraction(0) if abs(v) < snap else Fraction(v).limit_denominator(max_den) for v in res.x]
        ok, L = _exact_check(Phi, mu_counts, z, y)
        if ok:
            return {"pattern": int(z), "dual_value": float(res.fun), "verified": True, "denominator": int(L),
                    "nonzero_coefficients": int(sum(1 for q in y if q != 0))}
    return {"pattern": int(z), "dual_value": float(res.fun), "verified": False, "denominator": None}


def facial_set_from_support(Phi: np.ndarray, support: np.ndarray) -> np.ndarray:
    """Patterns positive under some distribution with the marginals of a distribution on `support`.

    One LP: maximize sum t_w over unobserved w subject to Phi^T d = 0, d_w >= t_w and
    0 <= t_w <= 1, with d free on observed patterns. A pattern w outside the support lies
    in the facial set exactly when some direction d that keeps all pairwise marginals
    fixed is positive at w and nonnegative off the support, so the optimum sets t_w = 1
    on the facial set and t_w = 0 elsewhere.
    """
    support = np.asarray(support, dtype=bool)
    K, m = Phi.shape
    U = np.flatnonzero(~support)
    out = support.copy()
    if len(U) == 0:
        return out
    nu = len(U)
    c = np.concatenate([np.zeros(K), -np.ones(nu)])
    A_eq = np.hstack([Phi.T.astype(float), np.zeros((m, nu))])
    A_ub = np.zeros((nu, K + nu))
    A_ub[np.arange(nu), U] = -1.0
    A_ub[np.arange(nu), K + np.arange(nu)] = 1.0
    bounds = [(None, None)] * K + [(0.0, 1.0)] * nu
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(nu), A_eq=A_eq, b_eq=np.zeros(m), bounds=bounds, method="highs")
    if res.status != 0:
        raise RuntimeError(f"facial-set LP failed: {res.message}")
    t = res.x[K:]
    out[U[t > 0.5]] = True
    return out


def analyse_counts(counts: np.ndarray, n: int, certify: bool = True) -> dict:
    """Pairwise hull, LP bounds and certificates for one table of pattern counts."""
    counts = np.asarray(counts, dtype=np.int64)
    N = int(counts.sum())
    Phi = moment_matrix(n)
    mu_counts = Phi.T @ counts
    mu = mu_counts / N
    S = {int(p) for p in np.flatnonzero(counts)}
    H = pairwise_hull(S, n)
    inside = sorted(H - S)
    bounds = {z: max_probability(Phi, mu, z) for z in inside}
    excluded = [z for z in inside if bounds[z] <= ZERO_TOL]
    possible = [z for z in inside if bounds[z] > ZERO_TOL]
    certs = [certificate(Phi, mu_counts, z) for z in excluded] if certify else []
    facial = facial_set_from_support(Phi, counts > 0)
    F_lp = set(S) | set(possible)
    pos_vals = [bounds[z] for z in possible]
    return {
        "issues": n, "respondents": N, "observed_patterns": len(S), "hull_patterns": len(H),
        "hull_unobserved": len(inside), "lp_possible_unobserved": len(possible),
        "lp_excluded_in_hull": len(excluded), "facial_set_size": len(F_lp),
        "certificates_verified": int(sum(c["verified"] for c in certs)),
        "certificates_attempted": len(certs),
        "max_denominator": max([c["denominator"] for c in certs if c["verified"]], default=None),
        "facial_set_agrees_with_support_lp": bool(F_lp == {int(p) for p in np.flatnonzero(facial)}),
        "smallest_positive_bound": min(pos_vals) if pos_vals else None,
        "largest_positive_bound": max(pos_vals) if pos_vals else None,
        "median_positive_bound": float(np.median(pos_vals)) if pos_vals else None,
        "hull_is_majority_closure": bool(H == majority_closure(S)),
        "excluded_patterns": excluded,
    }


# --- exhaustive comparison on small numbers of issues ---------------------------------

def hyperoctahedral_images(n: int) -> np.ndarray:
    """Images of every pattern under the 2^n n! permutations and flips of the items."""
    imgs = []
    for sigma in itertools.permutations(range(n)):
        for flip in range(1 << n):
            img = []
            for p in range(1 << n):
                q = 0
                for i in range(n):
                    q |= (((p >> i) & 1) ^ ((flip >> i) & 1)) << sigma[i]
                img.append(q)
            imgs.append(img)
    return np.array(imgs, dtype=np.int64)


def support_orbits(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Canonical representatives of nonempty supports under item permutations and flips."""
    imgs = hyperoctahedral_images(n)
    K = 1 << n
    masks = np.arange(1, 1 << K, dtype=np.int64)
    canon = np.full(len(masks), np.iinfo(np.int64).max, dtype=np.int64)
    for img in imgs:
        out = np.zeros_like(masks)
        for p in range(K):
            out |= ((masks >> p) & 1) << int(img[p])
        canon = np.minimum(canon, out)
    reps, sizes = np.unique(canon, return_counts=True)
    return reps, sizes


def _parity_rule_n3(S: set[int]) -> set[int]:
    """Closed form for three items: the facial set is S if the unobserved patterns have
    both parities, and all eight patterns otherwise (the kernel of the pairwise moment
    map is spanned by the parity character)."""
    U = [p for p in range(8) if p not in S]
    parities = {bin(p).count("1") % 2 for p in U}
    return set(S) if len(parities) == 2 else set(range(8))


def exhaustive_comparison(seed: int = 20261118) -> dict:
    """Compare LP exclusions with the pairwise hull for every support on three and four items.

    For each support, the moments come from integer weights drawn at random on the
    support; a second, independent weighting checks that the excluded set depends only
    on the support. Three items: all 255 supports, with the closed form as a third method.
    Four items: all 65,535 supports through their 401 orbits under permutations and flips
    of the items, which preserve the pairwise hull and the facial set.
    """
    rng = random.Random(seed)
    out = {}
    for n in (3, 4):
        K = 1 << n
        Phi = moment_matrix(n)
        if n == 3:
            reps = np.arange(1, 1 << K, dtype=np.int64)
            sizes = np.ones(len(reps), dtype=np.int64)
        else:
            reps, sizes = support_orbits(n)
        tot = int(sizes.sum())
        agree_methods = agree_weights = agree_closed = closed_checked = 0
        differs = strict_hull = 0
        unobserved_in_hull = excluded_in_hull = 0
        max_excluded = 0
        certs_ok = certs_total = 0
        example = None
        for mask, w in zip(reps, sizes):
            S = {p for p in range(K) if (int(mask) >> p) & 1}
            H = pairwise_hull(S, n)
            res = []
            for _ in range(2):
                counts = np.zeros(K, dtype=np.int64)
                for p in S:
                    counts[p] = rng.randint(1, 9)
                mu_counts = Phi.T @ counts
                mu = mu_counts / counts.sum()
                inside = sorted(H - S)
                bnd = {z: max_probability(Phi, mu, z) for z in inside}
                F = set(S) | {z for z in inside if bnd[z] > ZERO_TOL}
                res.append((F, counts, mu_counts, bnd))
            F, counts, mu_counts, bnd = res[0]
            agree_weights += int(w) * int(F == res[1][0])
            facial = facial_set_from_support(Phi, counts > 0)
            agree_methods += int(w) * int(F == {int(p) for p in np.flatnonzero(facial)})
            if n == 3:
                closed_checked += 1
                agree_closed += int(F == _parity_rule_n3(S))
            excl = sorted(H - F)
            for z in excl:
                c = certificate(Phi, mu_counts, z)
                certs_total += int(w)
                certs_ok += int(w) * int(c["verified"])
            strict_hull += int(w) * int(H != S)
            differs += int(w) * int(H != F)
            unobserved_in_hull += int(w) * len(H - S)
            excluded_in_hull += int(w) * len(excl)
            max_excluded = max(max_excluded, len(excl))
            if n == 3 and example is None and sorted(S) == [1, 2, 3, 4, 5, 6]:
                example = {"support": sorted(S), "hull": sorted(H), "facial_set": sorted(F),
                           "excluded": excl,
                           "bounds": {str(z): bnd[z] for z in sorted(bnd)}}
        out[f"issues_{n}"] = {
            "supports": tot, "methods_compared_on": len(reps),
            "facial_set_lp_agrees_per_pattern_lp": agree_methods,
            "same_exclusions_under_second_weighting": agree_weights,
            "closed_form_agrees": agree_closed if n == 3 else None,
            "closed_form_checked": closed_checked if n == 3 else None,
            "supports_with_unobserved_patterns_in_hull": strict_hull,
            "supports_where_lp_excludes_more_than_hull": differs,
            "unobserved_patterns_in_hull": unobserved_in_hull,
            "lp_excluded_patterns_in_hull": excluded_in_hull,
            "max_lp_excluded_in_hull_for_one_support": max_excluded,
            "certificates_verified": certs_ok, "certificates_total": certs_total,
        }
        if example is not None:
            out["six_profile_example"] = example
    out["all_agree"] = all(
        v["facial_set_lp_agrees_per_pattern_lp"] == v["supports"]
        and v["same_exclusions_under_second_weighting"] == v["supports"]
        and v["certificates_verified"] == v["certificates_total"]
        and (v["closed_form_agrees"] is None or v["closed_form_agrees"] == v["closed_form_checked"])
        for k, v in out.items() if k.startswith("issues_"))
    return out
