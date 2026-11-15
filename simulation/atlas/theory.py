"""Exhaustive verification of obstruction orders for representative commitment classes."""
from __future__ import annotations

import itertools
import random
import time

from . import classes as C
from . import constructions as K
from .geometry import (BoxClosure, HullClosure, convex_hull, largest_empty_convex_polygon,
                       random_general_position)
from .setsys import (Closure, bundled, full_mask, is_minimal_incompatible, order_by_families,
                     order_by_points)

SEED_PLANAR = 20261115


def boolean_orders(n_max: int = 5, n_max_unrestricted: int = 4, n_max_families: int = 3) -> dict:
    rows = []
    for name, gen in C.BOOLEAN_CLASSES.items():
        for n in range(1, n_max + 1):
            if name == "unrestricted" and n > n_max_unrestricted:
                continue
            N = 1 << n
            full = full_mask(N)
            atoms = gen(n)
            t0 = time.time()
            res = order_by_points(N, Closure(atoms, full))
            row = {"class": name, "issues": n, "atoms": len(set(atoms)), "order_by_points": res.order,
                   "nodes": res.nodes, "seconds": round(time.time() - t0, 3)}
            if n >= 2:
                row["formula"] = C.boolean_formula(name, n)
                row["matches_formula"] = res.order == row["formula"]
            if n <= n_max_families:
                row["order_by_families"] = order_by_families(atoms, full)
                row["methods_agree"] = row["order_by_families"] == res.order
            rows.append(row)
    return {"rows": rows,
            "all_formulas_match": all(r.get("matches_formula", True) for r in rows),
            "all_methods_agree": all(r.get("methods_agree", True) for r in rows),
            "rows_checked_against_formula": sum(1 for r in rows if "formula" in r),
            "rows_cross_checked": sum(1 for r in rows if "order_by_families" in r)}


def bundling_invariance(n: int = 2) -> dict:
    """Order of conjunction-closed classes equals the order of their atoms (by family enumeration)."""
    N = 1 << n
    full = full_mask(N)
    rows = []
    for name, gen in C.BOOLEAN_CLASSES.items():
        atoms = gen(n)
        bundles = bundled(atoms, full, max_arity=len(set(atoms)))
        rows.append({"class": name, "atoms": len(set(atoms)), "bundled_sets": len(bundles),
                     "order_atoms": order_by_families(atoms, full),
                     "order_bundled": order_by_families(bundles, full)})
    return {"issues": n, "rows": rows, "classes": len(rows),
            "invariant": all(r["order_atoms"] == r["order_bundled"] for r in rows)}


def witnesses(n_max: int = 10) -> dict:
    out = {"implication_chain": [], "parity": [], "double_cycle": [], "exclusion_triple": [],
           "unrestricted": [], "nominal_exclusions": [], "grid_cube": [], "arc_consistency_chain": []}
    for n in range(2, n_max + 1):
        full = full_mask(1 << n)
        out["implication_chain"].append({"issues": n, "size": n + 1,
                                         "minimal_incompatible": is_minimal_incompatible(K.implication_chain(n), full)})
        out["parity"].append({"issues": n, "size": n + 1,
                              "minimal_incompatible": is_minimal_incompatible(K.parity_family(n), full)})
        out["double_cycle"].append({"issues": n, "size": 2 * n,
                                    "minimal_incompatible": is_minimal_incompatible(K.double_cycle(n), full)})
        out["exclusion_triple"].append({"issues": n, "size": 3,
                                        "minimal_incompatible": is_minimal_incompatible(K.exclusion_triple(n), full)})
        out["arc_consistency_chain"].append(K.arc_consistency_detects_chain(n))
        if n <= 6:
            out["unrestricted"].append({"issues": n, "size": 2 ** n,
                                        "minimal_incompatible": is_minimal_incompatible(K.unrestricted_family(n), full)})
    for k in range(2, 11):
        out["nominal_exclusions"].append({"options": k, "size": k,
                                          "minimal_incompatible": is_minimal_incompatible(K.nominal_exclusion_family(k), full_mask(k))})
    for d in range(1, 5):
        out["grid_cube"].append(K.cube_family_on_grid(d, levels=3))
    ok = all(r.get("minimal_incompatible", True) for key in out for r in out[key] if isinstance(r, dict))
    ok &= all(r["vertices_separated"] for r in out["grid_cube"])
    ok &= all(r["domain_wipeout"] for r in out["arc_consistency_chain"])
    out["all_verified"] = bool(ok)
    return out


def parity_pairs(max_commitments: int = 12) -> dict:
    rows = [K.matched_parity_pair(c) for c in range(3, max_commitments + 1)]
    return {"rows": rows,
            "total_proper_subfamily_comparisons": sum(r["proper_subfamily_comparisons"] for r in rows),
            "all_equal": all(r["equal_counts_all_proper_subfamilies"] for r in rows),
            "all_even_compatible": all(r["even_compatible"] for r in rows),
            "all_odd_minimal": all(r["odd_minimal_incompatible"] for r in rows),
            "pairwise_jaccard_values": sorted({v for r in rows for v in r["pairwise_jaccard_values"]})}


def single_issue_and_products() -> dict:
    rows = []
    for k in range(2, 8):
        full = full_mask(k)
        rows.append({"case": "nominal_exclusions", "options": k,
                     "order": order_by_points(k, Closure(C.nominal_exclusions(k), full)).order, "expected": k})
        rows.append({"case": "nominal_choices", "options": k,
                     "order": order_by_points(k, Closure(C.nominal_choices(k), full)).order, "expected": 2})
        rows.append({"case": "ordinal_intervals", "options": k,
                     "order": order_by_points(k, Closure(C.ordinal_intervals(k), full)).order, "expected": 2})
    for sizes in ([3, 3], [2, 3, 4], [5, 5]):
        pts = C.product_space(sizes)
        full = full_mask(len(pts))
        excl = C.product_class(sizes, [C.nominal_exclusions(s) for s in sizes])
        intv = C.product_class(sizes, [C.ordinal_intervals(s) for s in sizes])
        rows.append({"case": "product_exclusions", "sizes": sizes,
                     "order": order_by_points(len(pts), Closure(excl, full)).order, "expected": max(sizes)})
        rows.append({"case": "product_intervals", "sizes": sizes,
                     "order": order_by_points(len(pts), Closure(intv, full)).order, "expected": 2})
    for L in range(2, 6):
        grid = list(itertools.product(range(L), repeat=2))
        rows.append({"case": "grid_boxes", "levels": L,
                     "order": order_by_points(len(grid), BoxClosure(grid)).order, "expected": 2})
        rows.append({"case": "grid_convex", "levels": L,
                     "order": order_by_points(len(grid), HullClosure(grid)).order, "expected": 4})
    return {"rows": rows, "cases": len(rows), "all_match": all(r["order"] == r["expected"] for r in rows)}


def planar_maps(m_values=range(4, 13), sets_per_m: int = 60, seed: int = SEED_PLANAR) -> dict:
    """Convex commitments over m occupied positions in general position in the plane."""
    rng = random.Random(seed)
    rows = []
    agree = 0
    total = 0
    for m in m_values:
        orders = []
        for _ in range(sets_per_m):
            pts = random_general_position(m, rng)
            h = order_by_points(m, HullClosure(pts)).order
            e, _ = largest_empty_convex_polygon(pts)
            orders.append(h)
            agree += int(h == e)
            total += 1
        hist = {str(v): orders.count(v) for v in sorted(set(orders))}
        rows.append({"positions": m, "sets": sets_per_m, "min_order": min(orders), "max_order": max(orders),
                     "mean_order": sum(orders) / len(orders), "histogram": hist})
    convex_rows = []
    for m in range(4, 11):
        # integer points on the parabola y = x^2 are in strictly convex position
        pts = [(x, x * x) for x in range(m)]
        assert len(convex_hull(pts)) == m
        convex_rows.append({"positions": m, "order": order_by_points(m, HullClosure(pts)).order})
    return {"seed": seed, "rows": rows, "sets_total": total,
            "order_equals_largest_empty_polygon": agree,
            "convex_position": convex_rows,
            "convex_position_order_equals_m": all(r["order"] == r["positions"] for r in convex_rows)}


def run(planar_m=range(4, 13), planar_sets: int = 60) -> dict:
    t0 = time.time()
    out = {
        "boolean": boolean_orders(),
        "bundling": bundling_invariance(2),
        "witnesses": witnesses(),
        "parity_pairs": parity_pairs(),
        "single_issue_and_products": single_issue_and_products(),
        "planar": planar_maps(m_values=planar_m, sets_per_m=planar_sets),
        "background_logic": background_logic_checks(),
        "pairwise_hull": pairwise_hull_checks(),
    }
    out["seconds"] = round(time.time() - t0, 1)
    return out


# --- background logic and the pairwise hull -----------------------------------------

def _majority(a: int, b: int, c: int) -> int:
    return (a & b) | (a & c) | (b & c)


def majority_closure(support: set[int]) -> set[int]:
    out = set(support)
    changed = True
    while changed:
        changed = False
        cur = list(out)
        for a, b, c in itertools.combinations_with_replacement(cur, 3):
            m = _majority(a, b, c)
            if m not in out:
                out.add(m)
                changed = True
    return out


def is_majority_closed(support: set[int]) -> bool:
    return all(_majority(a, b, c) in support for a, b, c in itertools.product(support, repeat=3))


def pairwise_hull(support: set[int], n: int) -> set[int]:
    """Positions whose every pair of answers occurs together in the support."""
    proj = {(i, j): {((p >> i) & 1, (p >> j) & 1) for p in support} for i, j in itertools.combinations(range(n), 2)}
    return {p for p in range(1 << n)
            if all(((p >> i) & 1, (p >> j) & 1) in proj[(i, j)] for i, j in proj)}


def categorical_order_within(support: set[int], n: int) -> int:
    """Obstruction order of categorical commitments when only positions in X are admissible."""
    pts = sorted(support)
    index = {p: k for k, p in enumerate(pts)}
    full = full_mask(len(pts))
    atoms = []
    for i in range(n):
        for a in (0, 1):
            m = 0
            for p in pts:
                if ((p >> i) & 1) == a:
                    m |= 1 << index[p]
            atoms.append(m)
    return order_by_points(len(pts), Closure(atoms, full)).order


def background_logic_checks(seed: int = 20261116, sample_n4: int = 3000) -> dict:
    """Categorical commitments relative to an admissible set X are decided pairwise
    (order <= 2) exactly when X is majority-closed; checked for every nonempty X over
    3 issues and a random sample over 4 issues."""
    out = {}
    rng = random.Random(seed)
    for n, sets in ((3, None), (4, sample_n4)):
        N = 1 << n
        if sets is None:
            masks = range(1, 1 << N)
        else:
            masks = [rng.randrange(1, 1 << N) for _ in range(sets)]
        agree = total = closed = 0
        for mk in masks:
            X = {p for p in range(N) if (mk >> p) & 1}
            pairwise_ok = categorical_order_within(X, n) <= 2
            mc = is_majority_closed(X)
            agree += int(pairwise_ok == mc)
            closed += int(mc)
            total += 1
        out[f"issues_{n}"] = {"admissible_sets": total, "majority_closed": closed, "agreements": agree,
                              "exhaustive": sets is None}
    out["all_agree"] = all(v["agreements"] == v["admissible_sets"] for v in out.values() if isinstance(v, dict))
    return out


def pairwise_hull_checks(seed: int = 20261117, samples: int = 3000) -> dict:
    """The answer patterns allowed by the pairwise tables of a support equal its majority closure."""
    rng = random.Random(seed)
    res = {}
    for n in (3, 4, 5):
        N = 1 << n
        if n == 3:
            masks = range(1, 1 << N)
        else:
            masks = [rng.randrange(1, 1 << N) for _ in range(samples)]
        agree = total = strict = 0
        for mk in masks:
            S = {p for p in range(N) if (mk >> p) & 1}
            h = pairwise_hull(S, n)
            agree += int(h == majority_closure(S))
            strict += int(h != S)
            total += 1
        res[f"issues_{n}"] = {"supports": total, "agreements": agree,
                              "supports_with_unobserved_patterns_in_hull": strict, "exhaustive": n == 3}
    res["all_agree"] = all(v["agreements"] == v["supports"] for v in res.values() if isinstance(v, dict))
    return res
