"""Exact planar convex geometry on integer points (no floating point)."""
from __future__ import annotations

import itertools
import random
from typing import Sequence

Point = tuple[int, int]


def cross(o: Point, a: Point, b: Point) -> int:
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_hull(points: Sequence[Point]) -> list[Point]:
    """Strict convex hull vertices in counter-clockwise order (monotone chain)."""
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts
    lower: list[Point] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: list[Point] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def in_hull(q: Point, hull: Sequence[Point]) -> bool:
    """Membership in the closed convex hull given by its ccw vertex list."""
    if not hull:
        return False
    if len(hull) == 1:
        return q == hull[0]
    if len(hull) == 2:
        a, b = hull
        if cross(a, b, q) != 0:
            return False
        return (min(a[0], b[0]) <= q[0] <= max(a[0], b[0])
                and min(a[1], b[1]) <= q[1] <= max(a[1], b[1]))
    for i in range(len(hull)):
        if cross(hull[i], hull[(i + 1) % len(hull)], q) < 0:
            return False
    return True


class HullClosure:
    """cl(A) = conv(A) intersected with a fixed finite point set S (S-convex hull)."""

    def __init__(self, points: Sequence[Point]):
        self.points = list(points)
        self._cache: dict[int, int] = {}

    def __call__(self, mask: int) -> int:
        hit = self._cache.get(mask)
        if hit is not None:
            return hit
        sub = [self.points[i] for i in range(len(self.points)) if (mask >> i) & 1]
        hull = convex_hull(sub)
        out = 0
        for i, q in enumerate(self.points):
            if in_hull(q, hull):
                out |= 1 << i
        self._cache[mask] = out
        return out


class BoxClosure:
    """cl(A) = smallest axis-parallel box containing A, intersected with S."""

    def __init__(self, points: Sequence[Point]):
        self.points = list(points)

    def __call__(self, mask: int) -> int:
        sub = [self.points[i] for i in range(len(self.points)) if (mask >> i) & 1]
        if not sub:
            return 0
        lo = [min(p[k] for p in sub) for k in range(2)]
        hi = [max(p[k] for p in sub) for k in range(2)]
        out = 0
        for i, q in enumerate(self.points):
            if all(lo[k] <= q[k] <= hi[k] for k in range(2)):
                out |= 1 << i
        return out


def largest_empty_convex_polygon(points: Sequence[Point]) -> tuple[int, list[Point]]:
    """Largest subset in strictly convex position whose closed hull holds no other point.

    Brute force over subsets in decreasing size; exact, for small point sets.
    """
    n = len(points)
    for r in range(n, 1, -1):
        for sub in itertools.combinations(range(n), r):
            pts = [points[i] for i in sub]
            hull = convex_hull(pts)
            if len(hull) != r:
                continue
            chosen = set(sub)
            if all(not in_hull(points[i], hull) for i in range(n) if i not in chosen):
                return r, pts
    return min(n, 1), list(points[:1])


def in_general_position(points: Sequence[Point]) -> bool:
    if len(set(points)) != len(points):
        return False
    for a, b, c in itertools.combinations(points, 3):
        if cross(a, b, c) == 0:
            return False
    return True


def random_general_position(m: int, rng: random.Random, scale: int = 10**6) -> list[Point]:
    while True:
        pts = [(rng.randrange(scale), rng.randrange(scale)) for _ in range(m)]
        if in_general_position(pts):
            return pts
