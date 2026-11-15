import hashlib
import itertools
import json
import random
from pathlib import Path

import numpy as np
import pytest

from atlas import classes as C
from atlas import constructions as K
from atlas import theory
from atlas.geometry import BoxClosure, HullClosure, convex_hull, in_hull, largest_empty_convex_polygon
from atlas.maxent import (MaxEnt, connected_information_3, design, g2_three_way, ipf_pairwise,
                          monomials, three_way_log_interaction)
from atlas.setsys import (Closure, full_mask, intersect, is_minimal_incompatible, jaccard_distance,
                          order_by_families, order_by_points)

DATA = Path(__file__).resolve().parents[1] / "data"


@pytest.mark.parametrize("name,n,expected", [
    ("categorical", 3, 2), ("implications_only", 3, 1), ("categorical_exclusions", 3, 3),
    ("categorical_implications", 3, 4), ("horn", 3, 4), ("parity", 3, 4),
    ("pairwise_clauses", 3, 6), ("unrestricted", 3, 8)])
def test_boolean_orders_small(name, n, expected):
    N = 1 << n
    full = full_mask(N)
    atoms = C.BOOLEAN_CLASSES[name](n)
    assert order_by_points(N, Closure(atoms, full)).order == expected
    assert order_by_families(atoms, full) == expected


def test_points_and_families_agree_on_random_set_systems():
    rng = random.Random(3)
    for _ in range(300):
        n_points = rng.randint(3, 7)
        full = full_mask(n_points)
        atoms = [rng.randrange(1, 1 << n_points) for _ in range(rng.randint(2, 8))]
        assert order_by_points(n_points, Closure(atoms, full)).order == order_by_families(atoms, full)


def test_witness_is_minimal_and_attains_order():
    res = order_by_points(1 << 3, Closure(C.two_clauses(3), full_mask(8)))
    assert is_minimal_incompatible(res.witness_family, full_mask(8))
    assert len(res.witness_family) == res.order == 6


@pytest.mark.parametrize("n", range(2, 8))
def test_constructions(n):
    full = full_mask(1 << n)
    assert is_minimal_incompatible(K.implication_chain(n), full)
    assert is_minimal_incompatible(K.parity_family(n), full)
    assert is_minimal_incompatible(K.double_cycle(n), full)
    assert len(K.double_cycle(n)) == 2 * n
    assert is_minimal_incompatible(K.exclusion_triple(n), full)


def test_grid_cube_family():
    for d in (1, 2, 3):
        r = K.cube_family_on_grid(d, 3)
        assert r["minimal_incompatible"] and r["vertices_separated"] and r["family_size"] == 2 ** d


def test_parity_pair_matches_pairwise_and_subfamily_counts():
    r = K.matched_parity_pair(5)
    assert r["equal_counts_all_proper_subfamilies"]
    assert r["pairwise_jaccard_values"] == ["2/3"]
    assert r["even_compatible"] and r["odd_minimal_incompatible"]


def test_jaccard_detects_pairwise_compatibility():
    assert jaccard_distance(0b0011, 0b1100) == 1.0
    assert jaccard_distance(0b0011, 0b0110) < 1.0


def test_geometry_hull_and_membership():
    pts = [(0, 0), (4, 0), (4, 4), (0, 4), (2, 2)]
    hull = convex_hull(pts)
    assert len(hull) == 4
    assert in_hull((2, 2), hull) and in_hull((4, 2), hull) and not in_hull((5, 2), hull)


def test_grid_orders():
    grid = list(itertools.product(range(3), repeat=2))
    assert order_by_points(len(grid), HullClosure(grid)).order == 4
    assert order_by_points(len(grid), BoxClosure(grid)).order == 2


def test_planar_order_equals_largest_empty_polygon():
    rng = random.Random(11)
    from atlas.geometry import random_general_position
    for _ in range(40):
        pts = random_general_position(rng.randint(4, 9), rng)
        assert order_by_points(len(pts), HullClosure(pts)).order == largest_empty_convex_polygon(pts)[0]


def test_arc_consistency_refutes_chain():
    for n in range(2, 9):
        assert K.arc_consistency_detects_chain(n)["domain_wipeout"]


def test_ipf_matches_pairwise_margins_and_xor():
    rng = np.random.default_rng(0)
    t = rng.integers(1, 50, size=(2, 2, 2)).astype(float)
    q = ipf_pairwise(t)
    p = t / t.sum()
    for ax in range(3):
        assert np.allclose(q.sum(axis=ax), p.sum(axis=ax), atol=1e-10)
    xor = np.zeros((2, 2, 2))
    for a, b, c in itertools.product((0, 1), repeat=3):
        if (a + b + c) % 2 == 0:
            xor[a, b, c] = 1.0
    ci = connected_information_3(xor)
    assert abs(ci["multi_information_bits"] - 1.0) < 1e-9
    assert abs(ci["third_order_bits"] - 1.0) < 1e-9
    assert abs(ci["pairwise_bits"]) < 1e-9


def test_g2_zero_for_pairwise_table():
    q = ipf_pairwise(np.array([[[10, 20], [30, 5]], [[7, 9], [11, 40]]], dtype=float)) * 1000
    g2, p = g2_three_way(q)
    assert g2 < 1e-6 and p > 0.99
    assert abs(three_way_log_interaction(q)) < 1e-6


def test_maxent_reproduces_margins():
    rng = np.random.default_rng(1)
    X = (rng.random((400, 4)) < 0.4).astype(float)
    X[:, 1] = np.where(rng.random(400) < 0.8, X[:, 0], X[:, 1])
    m = MaxEnt(4, 2).fit(X, lam=1e-8)
    p = m.probabilities()
    D = design(m.S, monomials(4, 2))
    assert np.allclose(p @ D, design(X, monomials(4, 2)).mean(axis=0), atol=1e-4)


def test_background_logic_and_hull_small():
    assert theory.is_majority_closed({0b000, 0b001, 0b011, 0b111})
    assert not theory.is_majority_closed({0b000, 0b011, 0b101, 0b111})
    assert not theory.is_majority_closed({0b011, 0b101, 0b110})
    S = {0b000, 0b011, 0b101, 0b110}
    assert theory.pairwise_hull(S, 3) == set(range(8)) == theory.majority_closure(S)


def test_frozen_data_match_manifest():
    manifest = json.loads((DATA / "MANIFEST.json").read_text())
    for name, rec in manifest["derived"].items():
        digest = hashlib.sha256((DATA / name).read_bytes()).hexdigest()
        assert digest == rec["sha256"], name


def test_exact_no_three_way_fit_matches_ipf_and_handles_boundary():
    from atlas.maxent import g2_batch, g2_fast, no3way_fit
    rng = np.random.default_rng(7)
    T = rng.multinomial(223, rng.dirichlet(np.ones(8)), size=500).reshape(-1, 2, 2, 2).astype(float)
    assert np.allclose(g2_fast(T), g2_batch(T), atol=1e-8)
    m = no3way_fit(T)
    assert np.allclose(m.sum(axis=1), T.sum(axis=1)) and np.allclose(m.sum(axis=(1, 2)), T.sum(axis=(1, 2)))
    inner = m[(m > 0).all(axis=(1, 2, 3))]
    assert len(inner) > 400
    lam = np.log(inner[:, 0, 0, 0] * inner[:, 0, 1, 1] * inner[:, 1, 0, 1] * inner[:, 1, 1, 0]) - np.log(
        inner[:, 0, 0, 1] * inner[:, 0, 1, 0] * inner[:, 1, 0, 0] * inner[:, 1, 1, 1])
    assert np.max(np.abs(lam)) < 1e-9
    # an empty cell of each parity: the table is its own fit
    t = np.array([9, 29, 3, 0, 0, 6, 9, 4], dtype=float).reshape(2, 2, 2)
    assert g2_fast(t[None])[0] == 0.0


def test_lp_certifies_exclusions_inside_the_hull():
    from atlas import lpbounds as L
    six = np.zeros(8, dtype=np.int64)
    six[[1, 2, 3, 4, 5, 6]] = 1
    r = L.analyse_counts(six, 3)
    assert r["hull_patterns"] == 8 and r["excluded_patterns"] == [0, 7]
    assert r["certificates_verified"] == 2 and r["facial_set_agrees_with_support_lp"]
    xor = np.zeros(8, dtype=np.int64)
    xor[[0, 3, 5, 6]] = 1
    r = L.analyse_counts(xor, 3)
    assert r["lp_excluded_in_hull"] == 0 and r["lp_possible_unobserved"] == 4
    assert abs(r["largest_positive_bound"] - 0.25) < 1e-9


def test_lp_three_items_closed_form():
    from atlas import lpbounds as L
    Phi = L.moment_matrix(3)
    for mask in range(1, 256):
        S = {p for p in range(8) if (mask >> p) & 1}
        counts = np.array([1 if p in S else 0 for p in range(8)], dtype=np.int64)
        F = L.facial_set_from_support(Phi, counts > 0)
        assert {int(p) for p in np.flatnonzero(F)} == L._parity_rule_n3(S)


def test_corrected_cv_variance():
    from atlas.empirical import corrected_cv_summary
    d = np.array([1.0, 2.0, 3.0, 4.0])
    r = corrected_cv_summary(d, np.array([10, 10, 10, 10]), 100)
    assert abs(r["se_corrected"] - np.sqrt((1 / 4 + 10 / 90) * d.var(ddof=1))) < 1e-12


def test_planted_model_keeps_margins():
    from atlas.empirical import planted_distribution
    rng = np.random.default_rng(3)
    X = (rng.random((200, 5)) < 0.5).astype(float)
    base = MaxEnt(5, 2).fit(X, 0.01)
    p, err = planted_distribution(base, (0, 2, 4), 1.5)
    assert err < 1e-6 and abs(p.sum() - 1) < 1e-12


def test_survey_screening_reproduces_source_counts():
    from atlas.empirical import screening_counts
    c = screening_counts(DATA / "consciousness_items.csv")
    assert c["submissions"] == 307 and c["eligible"] == 232 and c["eligible_attended_assc"] == 166
