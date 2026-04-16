#!/usr/bin/env python3
"""boundary_tri.py — boundary triangulation of a cannonball-lattice region.

Region: {(a,b,c) : 0<=a<=k, 0<=b<=l, 0<=c<=m, -r<=a+b-c<=q}
(5 parameters: k, l, m, r, q; all >=0, with r >= 0, q >= 1 typical)

Output: list of boundary triangles (as vertex tuples) tiling each of
the 8 bounding planes with a consistent template in the Cayley basis.

Includes BFS search for all parameter tuples with v <= VMAX, having
no degree-3 vertex and at least one face-interior degree-6 vertex.
"""
from __future__ import annotations

from collections import deque
from functools import lru_cache


Z = (0, 0, 0)
A = (1, 0, 0)
B = (0, 1, 0)
C = (0, 0, 1)
AC = (1, 0, 1)
BC = (0, 1, 1)
AMB = (1, -1, 0)


def add(p, q):
    return (p[0] + q[0], p[1] + q[1], p[2] + q[2])


def inside(x, k, l, m, r, q):
    a, b, c = x
    s = a + b - c
    return 0 <= a <= k and 0 <= b <= l and 0 <= c <= m and -r <= s <= q


def on_boundary(x, k, l, m, r, q):
    a, b, c = x
    s = a + b - c
    return (
        a == 0 or a == k or b == 0 or b == l or c == 0 or c == m
        or s == -r or s == q
    )


def cycle_to_smallest(tri):
    rots = [tri, (tri[1], tri[2], tri[0]), (tri[2], tri[0], tri[1])]
    return min(rots)


def reverse_orientation(tri):
    return (tri[0], tri[2], tri[1])


def facet_templates_lower():
    return {
        "a=0": [(Z, B, BC), (Z, BC, C)],
        "b=0": [(Z, C, AC), (Z, AC, A)],
        "c=0": [(Z, A, B), (Z, AMB, A)],
        "s=-r": [(Z, BC, AC), (Z, AC, AMB)],
    }


def facet_specs(k, l, m, r, q):
    lower = facet_templates_lower()
    def rev(name):
        return [reverse_orientation(t) for t in lower[name]]
    return [
        ("a=0", lambda x: x[0] == 0, lower["a=0"]),
        ("a=k", lambda x: x[0] == k, rev("a=0")),
        ("b=0", lambda x: x[1] == 0, lower["b=0"]),
        ("b=l", lambda x: x[1] == l, rev("b=0")),
        ("c=0", lambda x: x[2] == 0, lower["c=0"]),
        ("c=m", lambda x: x[2] == m, rev("c=0")),
        ("s=-r", lambda x: x[0] + x[1] - x[2] == -r, lower["s=-r"]),
        ("s=q", lambda x: x[0] + x[1] - x[2] == q, rev("s=-r")),
    ]


@lru_cache(maxsize=None)
def boundary_vertices(k, l, m, r, q):
    out = []
    for a in range(k + 1):
        for b in range(l + 1):
            for c in range(m + 1):
                x = (a, b, c)
                if inside(x, k, l, m, r, q) and on_boundary(x, k, l, m, r, q):
                    out.append(x)
    return tuple(out)


@lru_cache(maxsize=None)
def boundary_triangulation(k, l, m, r, q):
    verts = boundary_vertices(k, l, m, r, q)
    vertset = set(verts)
    tris = []
    seen = set()
    for _, is_on_facet, templates in facet_specs(k, l, m, r, q):
        for x in verts:
            if not is_on_facet(x):
                continue
            for t in templates:
                tri = tuple(add(x, u) for u in t)
                if any(v not in vertset for v in tri):
                    continue
                if not all(inside(v, k, l, m, r, q) for v in tri):
                    continue
                if not all(is_on_facet(v) for v in tri):
                    continue
                key = tuple(sorted(tri))
                if key in seen:
                    continue
                seen.add(key)
                tris.append(cycle_to_smallest(tri))
    tris.sort()
    first = ((0, 0, 0), (1, 0, 0), (0, 1, 0))
    if first in tris:
        tris.remove(first)
        tris.insert(0, first)
    return tuple(tris)


@lru_cache(maxsize=None)
def vertex_set_of_triangulation(k, l, m, r, q):
    tris = boundary_triangulation(k, l, m, r, q)
    return tuple(sorted({v for tri in tris for v in tri}))


@lru_cache(maxsize=None)
def graph_degrees(k, l, m, r, q):
    tris = boundary_triangulation(k, l, m, r, q)
    nbrs = {}
    for tri in tris:
        u, v, w = tri
        for x, y in [(u, v), (v, w), (w, u)]:
            nbrs.setdefault(x, set()).add(y)
            nbrs.setdefault(y, set()).add(x)
    return tuple((v, len(nbrs[v])) for v in sorted(nbrs))


def active_constraint_count(x, k, l, m, r, q):
    a, b, c = x
    s = a + b - c
    return sum([a == 0, a == k, b == 0, b == l, c == 0, c == m, s == -r, s == q])


def has_face_interior_degree6(k, l, m, r, q):
    for v, d in graph_degrees(k, l, m, r, q):
        if d == 6 and active_constraint_count(v, k, l, m, r, q) == 1:
            return True
    return False


def has_no_degree3(k, l, m, r, q):
    return all(d != 3 for _, d in graph_degrees(k, l, m, r, q))


def canonical_state(k, l, m, r, q):
    if k < 1 or l < 1 or m < 1 or r < 0 or q < 1:
        return None
    while True:
        old = (k, l, m, r, q)
        if k > l: k, l = l, k
        q = min(q, k + l)
        r = min(r, m)
        m = min(m, k + l + r)
        k = min(k, q + m)
        l = min(l, q + m)
        if k > l: k, l = l, k
        if old == (k, l, m, r, q):
            return (k, l, m, r, q)


def neighbors(state):
    k, l, m, r, q = state
    raw = [
        (k + 1, l, m, r, q),
        (k, l + 1, m, r, q),
        (k, l, m + 1, r, q),
        (k, l, m, r + 1, q),
        (k, l, m, r, q + 1),
    ]
    out = []
    for s in raw:
        t = canonical_state(*s)
        if t is not None and t != state:
            out.append(t)
    return out


def search(vmax=30):
    fmax = 2 * vmax - 4
    start = canonical_state(1, 1, 1, 0, 1)
    q = deque([start])
    seen = {start}
    out = []
    while q:
        state = q.popleft()
        k, l, m, r, qq = state
        tris = boundary_triangulation(k, l, m, r, qq)
        f = len(tris)
        if f > fmax:
            continue
        verts = vertex_set_of_triangulation(k, l, m, r, qq)
        v = len(verts)
        if f == 2 * v - 4:
            if v <= vmax and has_no_degree3(k, l, m, r, qq) and has_face_interior_degree6(k, l, m, r, qq):
                out.append((state, v, f))
        for nxt in neighbors(state):
            if nxt not in seen:
                seen.add(nxt)
                q.append(nxt)
    out.sort(key=lambda x: (x[1], x[2], x[0]))
    return out


def main():
    import sys
    vmax = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    for state, v, f in search(vmax):
        print(state, "V=", v, "F=", f)


if __name__ == "__main__":
    main()
