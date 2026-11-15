"""Atomic commitment classes over finite position spaces.

Boolean classes live on the 2^n complete yes/no positions over n issues; point p
encodes the answer to issue i in bit i. Each function returns a list of atomic
commitments as bit masks over the 2^n points. Conjunctions of atoms (bundled
positions) never change the obstruction order, so only atoms are listed.
"""
from __future__ import annotations

import itertools
from typing import Callable

from .setsys import full_mask


def _mask(n_points: int, pred: Callable[[int], bool]) -> int:
    m = 0
    for p in range(n_points):
        if pred(p):
            m |= 1 << p
    return m


def _bit(p: int, i: int) -> int:
    return (p >> i) & 1


def literals(n: int) -> list[int]:
    """Categorical answers: issue i answered a."""
    N = 1 << n
    return [_mask(N, lambda p, i=i, a=a: _bit(p, i) == a) for i in range(n) for a in (0, 1)]


def implications(n: int) -> list[int]:
    """Conditional commitments x_i -> x_j between distinct issues."""
    N = 1 << n
    return [_mask(N, lambda p, i=i, j=j: (not _bit(p, i)) or _bit(p, j))
            for i in range(n) for j in range(n) if i != j]


def exclusions(n: int) -> list[int]:
    """Not both x_i and x_j (NAND)."""
    N = 1 << n
    return [_mask(N, lambda p, i=i, j=j: not (_bit(p, i) and _bit(p, j)))
            for i, j in itertools.combinations(range(n), 2)]


def two_clauses(n: int) -> list[int]:
    """All clauses with at most two literals (2-CNF atoms), units included."""
    N = 1 << n
    out = literals(n)
    for i, j in itertools.combinations(range(n), 2):
        for a in (0, 1):
            for b in (0, 1):
                out.append(_mask(N, lambda p, i=i, j=j, a=a, b=b: _bit(p, i) == a or _bit(p, j) == b))
    return out


def horn_clauses(n: int) -> list[int]:
    """Clauses with at most one positive literal (rules with conjunctive premises)."""
    N = 1 << n
    out = []
    for neg in range(1 << n):
        heads = [None] + [i for i in range(n) if not (neg >> i) & 1]
        for head in heads:
            if neg == 0 and head is None:
                continue
            out.append(_mask(N, lambda p, neg=neg, head=head:
                             any(_bit(p, i) == 0 for i in range(n) if (neg >> i) & 1)
                             or (head is not None and _bit(p, head) == 1)))
    return out


def parity_constraints(n: int) -> list[int]:
    """Affine hyperplanes a.x = b over GF(2), a nonzero."""
    N = 1 << n
    return [_mask(N, lambda p, a=a, b=b: (p & a).bit_count() % 2 == b)
            for a in range(1, 1 << n) for b in (0, 1)]


def unrestricted(n: int) -> list[int]:
    """Complements of single positions; their intersections give every subset."""
    N = 1 << n
    full = full_mask(N)
    return [full ^ (1 << p) for p in range(N)]


BOOLEAN_CLASSES: dict[str, Callable[[int], list[int]]] = {
    "categorical": literals,
    "implications_only": implications,
    "categorical_exclusions": lambda n: literals(n) + exclusions(n),
    "categorical_implications": lambda n: literals(n) + implications(n),
    "horn": horn_clauses,
    "parity": parity_constraints,
    "pairwise_clauses": two_clauses,
    "unrestricted": unrestricted,
}


def boolean_formula(name: str, n: int) -> int:
    """Proven obstruction order of each Boolean class on n >= 2 issues."""
    return {
        "categorical": 2,
        "implications_only": 1,
        "categorical_exclusions": 3,
        "categorical_implications": n + 1,
        "horn": n + 1,
        "parity": n + 1,
        "pairwise_clauses": 2 * n,
        "unrestricted": 2 ** n,
    }[name]


# --- one nominal or ordinal issue -------------------------------------------------

def nominal_exclusions(k: int) -> list[int]:
    """One issue with k exhaustive options; commitments 'not option a'."""
    full = full_mask(k)
    return [full ^ (1 << a) for a in range(k)]


def nominal_choices(k: int) -> list[int]:
    """One issue with k options; commitments 'option a'."""
    return [1 << a for a in range(k)]


def ordinal_intervals(levels: int) -> list[int]:
    """One ordinal issue (e.g. a Likert item); commitments are intervals [lo, hi]."""
    out = []
    for lo in range(levels):
        for hi in range(lo, levels):
            m = 0
            for v in range(lo, hi + 1):
                m |= 1 << v
            out.append(m)
    return out


def product_space(sizes: list[int]) -> list[tuple[int, ...]]:
    return list(itertools.product(*[range(s) for s in sizes]))


def product_class(sizes: list[int], per_issue: list[list[int]]) -> list[int]:
    """Commitments constraining one issue at a time (cylinders), on the product space."""
    pts = product_space(sizes)
    out = []
    for q, atoms in enumerate(per_issue):
        for a in atoms:
            m = 0
            for idx, p in enumerate(pts):
                if (a >> p[q]) & 1:
                    m |= 1 << idx
            out.append(m)
    return out
