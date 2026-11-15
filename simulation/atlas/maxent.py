"""Maximum-entropy models of binary position data and connected information.

For n binary items the k-th order maximum-entropy model is the distribution of
largest entropy that reproduces all marginals of order <= k. Its exponential-family
form has one parameter per monomial of degree <= k. The connected information of
order k is the entropy drop from the (k-1)-th to the k-th order model (Schneidman,
Still, Berry and Bialek 2003); the multi-information is their sum.

Three-item tables are fitted exactly by iterative proportional fitting of the three
two-way margins. Larger item sets are fitted by penalised maximum likelihood over
all 2^n states.
"""
from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

LOG2 = math.log(2.0)


def states(n: int) -> np.ndarray:
    """All 2^n binary vectors; column i is item i, row index = sum x_i 2^(n-1-i)."""
    return np.array(list(itertools.product((0, 1), repeat=n)), dtype=np.float64)


def monomials(n: int, order: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for k in range(1, order + 1):
        out.extend(itertools.combinations(range(n), k))
    return out


def design(X: np.ndarray, terms: list[tuple[int, ...]]) -> np.ndarray:
    cols = []
    for t in terms:
        col = np.ones(X.shape[0])
        for i in t:
            col = col * X[:, i]
        cols.append(col)
    return np.stack(cols, axis=1)


def entropy_bits(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=np.float64).ravel()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


# --- three-item tables ------------------------------------------------------------

def ipf_pairwise(table: np.ndarray, iters: int = 2000, tol: float = 1e-13) -> np.ndarray:
    """Fit of the no-three-way-interaction model to (..., 2, 2, 2) tables.

    Works on a batch of tables (leading axes). Returns probabilities that match all
    three two-way margins of the normalised input.
    """
    t = np.asarray(table, dtype=np.float64)
    p = t / t.sum(axis=(-3, -2, -1), keepdims=True)
    q = np.full_like(p, 1.0 / 8.0)
    for _ in range(iters):
        q_old = q
        for axis in (-1, -2, -3):
            target = p.sum(axis=axis, keepdims=True)
            current = q.sum(axis=axis, keepdims=True)
            ratio = np.divide(target, current, out=np.zeros_like(target), where=current > 0)
            q = q * ratio
        if np.max(np.abs(q - q_old)) < tol:
            break
    return q


def independent_fit(table: np.ndarray) -> np.ndarray:
    t = np.asarray(table, dtype=np.float64)
    p = t / t.sum()
    m0 = p.sum(axis=(1, 2))
    m1 = p.sum(axis=(0, 2))
    m2 = p.sum(axis=(0, 1))
    return m0[:, None, None] * m1[None, :, None] * m2[None, None, :]


def connected_information_3(table: np.ndarray) -> dict:
    """Multi-information and its split into pairwise and third-order parts (bits)."""
    t = np.asarray(table, dtype=np.float64)
    p = t / t.sum()
    h1 = entropy_bits(independent_fit(t))
    h2 = entropy_bits(ipf_pairwise(t))
    h3 = entropy_bits(p)
    multi = h1 - h3
    i3 = max(h2 - h3, 0.0)
    return {"multi_information_bits": multi, "pairwise_bits": h1 - h2,
            "third_order_bits": i3,
            "third_order_fraction": (i3 / multi) if multi > 0 else 0.0}


def g2_three_way(table: np.ndarray) -> tuple[float, float]:
    """Likelihood-ratio statistic for no three-way interaction (1 df) and p-value."""
    t = np.asarray(table, dtype=np.float64)
    n = t.sum()
    q = ipf_pairwise(t) * n
    mask = t > 0
    g2 = float(2.0 * np.sum(t[mask] * np.log(t[mask] / q[mask])))
    g2 = max(g2, 0.0)
    return g2, float(chi2.sf(g2, 1))


def g2_batch(tables: np.ndarray) -> np.ndarray:
    t = np.asarray(tables, dtype=np.float64)
    n = t.sum(axis=(-3, -2, -1), keepdims=True)
    q = ipf_pairwise(t) * n
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(t > 0, t * np.log(np.where(t > 0, t, 1.0) / np.where(q > 0, q, 1.0)), 0.0)
    return np.maximum(2.0 * term.sum(axis=(-3, -2, -1)), 0.0)


EVEN_CELLS = [0, 3, 5, 6]  # 000, 011, 101, 110 in C order (index 4a + 2b + c)
ODD_CELLS = [1, 2, 4, 7]   # 001, 010, 100, 111


def no3way_fit(tables: np.ndarray, iters: int = 200) -> np.ndarray:
    """Exact fit of the no-three-way-interaction model to (..., 2, 2, 2) count tables.

    The tables with the two-way margins of t are t + delta * s, where s is +1 on cells
    of even parity and -1 on cells of odd parity. The fit is the member of this line
    with zero three-way log interaction, the root of the strictly increasing function
    sum_even log(t + delta) - sum_odd log(t - delta) on (-min_even t, min_odd t), found
    by Newton steps kept inside a shrinking bracket. When both parities contain an
    empty cell the line meets the nonnegative orthant only at t, which is then the fit
    (the extended maximum-likelihood estimate that iterative scaling converges to).
    """
    t = np.asarray(tables, dtype=np.float64)
    flat = t.reshape(-1, 8)
    te, to = flat[:, EVEN_CELLS], flat[:, ODD_CELLS]
    lo, hi = -te.min(axis=1), to.min(axis=1)
    point = hi - lo <= 0.0
    x = np.where((lo < 0.0) & (hi > 0.0), 0.0, 0.5 * (lo + hi))
    x[point] = 0.0
    act = ~point
    for _ in range(iters):
        if not act.any():
            break
        a = te[act] + x[act, None]
        b = to[act] - x[act, None]
        g = np.log(a).sum(axis=1) - np.log(b).sum(axis=1)
        gp = (1.0 / a).sum(axis=1) + (1.0 / b).sum(axis=1)
        la, ha, xa = lo[act], hi[act], x[act]
        ha = np.where(g > 0.0, xa, ha)
        la = np.where(g < 0.0, xa, la)
        xn = xa - g / gp
        bad = ~((xn > la) & (xn < ha))
        xn = np.where(bad, 0.5 * (la + ha), xn)
        done = (np.abs(xn - xa) <= 1e-13 * (1.0 + np.abs(xa))) | (g == 0.0)
        lo[act], hi[act], x[act] = la, ha, np.where(g == 0.0, xa, xn)
        idx = np.flatnonzero(act)
        act[idx[done]] = False
    m = flat.copy()
    m[:, EVEN_CELLS] += x[:, None]
    m[:, ODD_CELLS] -= x[:, None]
    return m.reshape(t.shape)


def g2_fast(tables: np.ndarray) -> np.ndarray:
    """Likelihood-ratio statistic for no three-way interaction, from the exact fit."""
    t = np.asarray(tables, dtype=np.float64)
    m = no3way_fit(t)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(t > 0, t * np.log(np.where(t > 0, t, 1.0) / np.where(m > 0, m, 1.0)), 0.0)
    return np.maximum(2.0 * term.sum(axis=(-3, -2, -1)), 0.0)


def three_way_log_interaction(table: np.ndarray, pseudo: float = 0.0) -> float:
    """log of the ratio of conditional odds ratios, log[OR_12|3=1 / OR_12|3=0]."""
    t = np.asarray(table, dtype=np.float64) + pseudo
    lt = np.log(t)
    total = 0.0
    for a, b, c in itertools.product((0, 1), repeat=3):
        sign = (-1) ** (3 - (a + b + c))
        total += sign * lt[a, b, c]
    return float(total)


def table3(X: np.ndarray, i: int, j: int, k: int, weights: np.ndarray | None = None) -> np.ndarray:
    idx = (X[:, i] * 4 + X[:, j] * 2 + X[:, k]).astype(int)
    w = np.ones(len(idx)) if weights is None else weights
    return np.bincount(idx, weights=w, minlength=8).reshape(2, 2, 2)


# --- n-item maximum-entropy models -------------------------------------------------

class MaxEnt:
    """Penalised maximum-likelihood fit of the order-k maximum-entropy model.

    The objective is -mean log-likelihood + (lam / 2) * ||theta||^2. With lam -> 0
    the fit reproduces the empirical marginals of order <= k (the maximum-entropy
    solution); lam > 0 shrinks interactions towards zero.
    """

    def __init__(self, n: int, order: int):
        self.n = n
        self.order = order
        self.terms = monomials(n, order)
        self.S = states(n)
        self.F = design(self.S, self.terms)
        self.theta = np.zeros(len(self.terms))

    def log_partition(self, theta: np.ndarray) -> tuple[float, np.ndarray]:
        e = self.F @ theta
        m = e.max()
        w = np.exp(e - m)
        z = w.sum()
        return float(np.log(z) + m), w / z

    def fit(self, X: np.ndarray, lam: float, weights: np.ndarray | None = None) -> "MaxEnt":
        D = design(X, self.terms)
        w = np.ones(len(X)) if weights is None else np.asarray(weights, dtype=np.float64)
        mu = (w[:, None] * D).sum(axis=0) / w.sum()

        def objective(th: np.ndarray) -> tuple[float, np.ndarray]:
            logz, p = self.log_partition(th)
            val = -(mu @ th - logz) + 0.5 * lam * float(th @ th)
            grad = -(mu - p @ self.F) + lam * th
            return val, grad

        res = minimize(objective, self.theta.copy(), jac=True, method="L-BFGS-B",
                       options={"maxiter": 10000, "gtol": 1e-9})
        self.theta = res.x
        return self

    def probabilities(self) -> np.ndarray:
        return self.log_partition(self.theta)[1]

    def mean_loglik_bits(self, X: np.ndarray) -> float:
        logz, _ = self.log_partition(self.theta)
        D = design(X, self.terms)
        return float(((D @ self.theta) - logz).mean() / LOG2)

    def sample(self, n_rows: int, rng: np.random.Generator) -> np.ndarray:
        p = self.probabilities()
        idx = rng.choice(len(p), size=n_rows, p=p)
        return self.S[idx]


def state_index(X: np.ndarray) -> np.ndarray:
    n = X.shape[1]
    weights = 2 ** np.arange(n - 1, -1, -1)
    return (X @ weights).astype(int)


def marginal3_from_full(p_full: np.ndarray, n: int, i: int, j: int, k: int) -> np.ndarray:
    S = states(n)
    idx = (S[:, i] * 4 + S[:, j] * 2 + S[:, k]).astype(int)
    return np.bincount(idx, weights=p_full, minlength=8).reshape(2, 2, 2)


def benjamini_hochberg(pvals: np.ndarray, q: float) -> np.ndarray:
    p = np.asarray(pvals)
    m = len(p)
    order = np.argsort(p)
    thresh = q * (np.arange(1, m + 1) / m)
    passed = p[order] <= thresh
    k = np.max(np.nonzero(passed)[0]) + 1 if passed.any() else 0
    out = np.zeros(m, dtype=bool)
    out[order[:k]] = True
    return out
