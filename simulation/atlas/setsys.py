"""Exact computations on finite set systems.

A position space is a finite set of N points, indexed 0..N-1. A commitment is a
subset of the space, stored as a Python integer bit mask. A family of commitments
is compatible when the intersection of its members is nonempty.

The obstruction order of a class of commitments is the largest size of a minimal
incompatible family drawn from the class (a family with empty intersection all of
whose proper subfamilies have nonempty intersection). It equals the Helly number of
the class. Two independent procedures compute it:

* ``order_by_points``: closes the class under intersection and searches for the
  largest set F of points with empty intersection of the closures cl(F minus p)
  ("Helly-independent" sets). The search uses the fact that subsets of such sets
  have the same property, so it can prune.
* ``order_by_families``: enumerates every subfamily of the atomic commitments and
  keeps the largest minimal incompatible one. It is exponential in the number of
  atoms and is used only on small cases as a cross-check.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Callable, Iterable, Sequence


def popcount(x: int) -> int:
    return x.bit_count()


def full_mask(n_points: int) -> int:
    return (1 << n_points) - 1


def intersect(family: Iterable[int], full: int) -> int:
    out = full
    for member in family:
        out &= member
        if out == 0:
            break
    return out


def is_minimal_incompatible(family: Sequence[int], full: int) -> bool:
    """True if the family has empty intersection and every proper subfamily meets.

    Compatibility is inherited by subfamilies, so it suffices to check the
    subfamilies obtained by deleting one member.
    """
    if intersect(family, full) != 0:
        return False
    for k in range(len(family)):
        if intersect((m for j, m in enumerate(family) if j != k), full) == 0:
            return False
    return True


class Closure:
    """Closure operator of the intersection closure of a list of atomic sets.

    cl(S) is the intersection of all atoms containing S (the whole space if none
    does). Results are cached by mask.
    """

    def __init__(self, atoms: Iterable[int], full: int):
        self.atoms = sorted(set(atoms))
        self.full = full
        self._cache: dict[int, int] = {}

    def __call__(self, s: int) -> int:
        hit = self._cache.get(s)
        if hit is not None:
            return hit
        out = self.full
        for a in self.atoms:
            if s & ~a == 0:
                out &= a
        self._cache[s] = out
        return out


@dataclass
class OrderResult:
    order: int
    witness_points: list[int]
    witness_family: list[int]
    nodes: int


def helly_independent(points: Sequence[int], cl: Callable[[int], int], full: int) -> bool:
    """Empty intersection of cl(F minus p) over p in F (for |F| >= 2)."""
    inter = full
    for k in range(len(points)):
        sub = 0
        for j, p in enumerate(points):
            if j != k:
                sub |= 1 << p
        inter &= cl(sub)
        if inter == 0:
            return True
    return inter == 0


def order_by_points(n_points: int, cl: Callable[[int], int]) -> OrderResult:
    """Largest Helly-independent point set; equals the obstruction order.

    Returns order 1 when no two points form a Helly-independent pair, which
    means every finite family of commitments from the class is compatible.
    """
    full = full_mask(n_points)
    best: list = [1, [], 0]
    nodes = 0

    def dfs(chosen: list[int], start: int) -> None:
        nonlocal nodes
        for p in range(start, n_points):
            cand = chosen + [p]
            if len(cand) >= 2 and not helly_independent(cand, cl, full):
                continue
            nodes += 1
            if len(cand) > best[0]:
                best[0] = len(cand)
                best[1] = list(cand)
            dfs(cand, p + 1)

    dfs([], 0)
    order, pts = best[0], best[1]
    family: list[int] = []
    if order >= 2:
        for k in range(len(pts)):
            sub = 0
            for j, p in enumerate(pts):
                if j != k:
                    sub |= 1 << p
            family.append(cl(sub))
        assert is_minimal_incompatible(family, full)
    return OrderResult(order=order, witness_points=pts, witness_family=family, nodes=nodes)


def order_by_families(atoms: Iterable[int], full: int) -> int:
    """Largest minimal incompatible subfamily of distinct atoms, by full enumeration."""
    atoms = sorted(set(atoms))
    m = len(atoms)
    best = 1
    for mask in range(1, 1 << m):
        r = mask.bit_count()
        if r <= best:
            continue
        fam = [atoms[i] for i in range(m) if (mask >> i) & 1]
        if is_minimal_incompatible(fam, full):
            best = r
    return best


def bundled(atoms: Sequence[int], full: int, max_arity: int) -> list[int]:
    """All nonempty intersections of at most ``max_arity`` atoms (conjunctions)."""
    atoms = sorted(set(atoms))
    out = set(atoms)
    for r in range(2, max_arity + 1):
        for combo in itertools.combinations(atoms, r):
            out.add(intersect(combo, full))
    out.discard(0)
    return sorted(out)


def jaccard_distance(a: int, b: int) -> float:
    union = popcount(a | b)
    if union == 0:
        return 0.0
    return 1.0 - popcount(a & b) / union
