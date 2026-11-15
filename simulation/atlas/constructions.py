"""Explicit minimal incompatible families that attain the lower bounds."""
from __future__ import annotations

import itertools
from fractions import Fraction

from .setsys import full_mask, intersect, is_minimal_incompatible, popcount


def _mask(N, pred):
    m = 0
    for p in range(N):
        if pred(p):
            m |= 1 << p
    return m


def _b(p, i):
    return (p >> i) & 1


def implication_chain(n: int) -> list[int]:
    """x_1, x_1 -> x_2, ..., x_{n-1} -> x_n, not x_n: categorical plus implications (and Horn)."""
    N = 1 << n
    fam = [_mask(N, lambda p: _b(p, 0) == 1)]
    for i in range(n - 1):
        fam.append(_mask(N, lambda p, i=i: (not _b(p, i)) or _b(p, i + 1)))
    fam.append(_mask(N, lambda p: _b(p, n - 1) == 0))
    return fam


def parity_family(n: int) -> list[int]:
    """x_i = 0 for every issue, and x_1 + ... + x_n = 1 (mod 2): n + 1 parity constraints."""
    N = 1 << n
    fam = [_mask(N, lambda p, i=i: _b(p, i) == 0) for i in range(n)]
    fam.append(_mask(N, lambda p: p.bit_count() % 2 == 1))
    return fam


def double_cycle(n: int) -> list[int]:
    """2n two-literal clauses: x_1 -> ... -> x_n -> not x_1 and not x_1 -> ... -> not x_n -> x_1."""
    N = 1 << n
    fam = []
    for i in range(n - 1):
        fam.append(_mask(N, lambda p, i=i: (not _b(p, i)) or _b(p, i + 1)))      # x_i -> x_{i+1}
    fam.append(_mask(N, lambda p: (not _b(p, n - 1)) or (not _b(p, 0))))        # x_n -> not x_1
    for i in range(n - 1):
        fam.append(_mask(N, lambda p, i=i: _b(p, i) or (not _b(p, i + 1))))      # not x_i -> not x_{i+1}
    fam.append(_mask(N, lambda p: _b(p, n - 1) or _b(p, 0)))                    # not x_n -> x_1
    return fam


def exclusion_triple(n: int) -> list[int]:
    N = 1 << n
    return [_mask(N, lambda p: _b(p, 0) == 1), _mask(N, lambda p: _b(p, 1) == 1),
            _mask(N, lambda p: not (_b(p, 0) and _b(p, 1)))]


def unrestricted_family(n: int) -> list[int]:
    N = 1 << n
    full = full_mask(N)
    return [full ^ (1 << p) for p in range(N)]


def nominal_exclusion_family(k: int) -> list[int]:
    full = full_mask(k)
    return [full ^ (1 << a) for a in range(k)]


def cube_family_on_grid(d: int, levels: int) -> dict:
    """Doignon's lower bound: the 2^d sets {0,1}^d minus one vertex, on the grid {0..L-1}^d.

    Each set is grid-convex: the only grid points of [0,1]^d are cube vertices, and the
    removed vertex v is separated from the others by f(x) = sum_i (x_i if v_i else 1 - x_i),
    which equals d at v and at most d - 1 at every other vertex.
    """
    grid = list(itertools.product(range(levels), repeat=d))
    index = {g: i for i, g in enumerate(grid)}
    cube = list(itertools.product((0, 1), repeat=d))
    in_unit_box = [g for g in grid if all(c in (0, 1) for c in g)]
    assert sorted(in_unit_box) == sorted(cube)
    full = full_mask(len(grid))
    fam = []
    separated = True
    for v in cube:
        others = [c for c in cube if c != v]
        f = lambda x: sum(x[i] if v[i] else 1 - x[i] for i in range(d))
        separated &= (f(v) == d) and max(f(c) for c in others) <= d - 1
        m = 0
        for c in others:
            m |= 1 << index[c]
        fam.append(m)
    return {"dimension": d, "levels": levels, "grid_points": len(grid), "family_size": len(fam),
            "vertices_separated": bool(separated),
            "minimal_incompatible": is_minimal_incompatible(fam, full)}


def matched_parity_pair(n_commitments: int) -> dict:
    """The standard parity pair: equal pairwise Jaccard distances and subfamily counts.

    Issues: n_commitments - 1 binary issues. Commitments: x_i = 0 for each issue, plus
    x_1 + ... + x_m = b with b = 0 (compatible family) or b = 1 (minimal obstruction).
    """
    m = n_commitments - 1
    N = 1 << m
    full = full_mask(N)
    base = [_mask(N, lambda p, i=i: _b(p, i) == 0) for i in range(m)]
    even = base + [_mask(N, lambda p: p.bit_count() % 2 == 0)]
    odd = base + [_mask(N, lambda p: p.bit_count() % 2 == 1)]
    comparisons = 0
    equal_counts = True
    for r in range(1, n_commitments):
        for sub in itertools.combinations(range(n_commitments), r):
            ce = popcount(intersect((even[i] for i in sub), full))
            co = popcount(intersect((odd[i] for i in sub), full))
            equal_counts &= (ce == co == 2 ** (m - r))
            comparisons += 1
    pair_distances = set()
    for i, j in itertools.combinations(range(n_commitments), 2):
        for fam in (even, odd):
            a, b = fam[i], fam[j]
            pair_distances.add(Fraction(popcount(a ^ b), popcount(a | b)))
    return {"commitments": n_commitments,
            "proper_subfamily_comparisons": comparisons,
            "equal_counts_all_proper_subfamilies": bool(equal_counts),
            "pairwise_jaccard_values": sorted(str(x) for x in pair_distances),
            "even_compatible": intersect(even, full) != 0,
            "odd_minimal_incompatible": is_minimal_incompatible(odd, full)}


def arc_consistency_detects_chain(n: int) -> dict:
    """Issue-local propagation (unary and binary constraints) on the implication chain.

    Domains start at {0,1}; the unary commitments x_1 = 1 and x_n = 0 prune domains;
    each implication x_i -> x_{i+1} removes 0 from D(x_{i+1}) when D(x_i) = {1}.
    The chain is refuted when some domain empties, after checking pairs of issues only.
    """
    dom = [{0, 1} for _ in range(n)]
    dom[0] &= {1}
    dom[n - 1] &= {0}
    changed = True
    steps = 0
    while changed:
        changed = False
        for i in range(n - 1):
            # support for x_{i+1} = v requires some u in D(x_i) with (u -> v)
            new_next = {v for v in dom[i + 1] if any((not u) or v for u in dom[i])}
            new_prev = {u for u in dom[i] if any((not u) or v for v in dom[i + 1])}
            if new_next != dom[i + 1] or new_prev != dom[i]:
                dom[i + 1], dom[i] = new_next, new_prev
                changed = True
                steps += 1
    wiped = any(len(d) == 0 for d in dom)
    return {"issues": n, "commitments": n + 1, "domain_wipeout": wiped,
            "largest_scope_examined": 2, "propagation_steps": steps}
