#!/usr/bin/env python3
"""boundary_tri.py — boundary triangulation of a cannonball-lattice region.

(From ChatGPT. Input: five integers k, l, m, p, q defining the octahedral
region {(a,b,c) : 0<=a<=k, 0<=b<=l, 0<=c<=m, p<=a+b-c<=q}.)

Output: list of boundary triangles (as vertex tuples) that tile each
of the 8 bounding faces with a consistent template.
"""
from __future__ import annotations

import sys


Z = (0, 0, 0)
A = (1, 0, 0)
B = (0, 1, 0)
C = (0, 0, 1)
AC = (1, 0, 1)
BC = (0, 1, 1)
AMB = (1, -1, 0)


def add(p, q):
    return (p[0] + q[0], p[1] + q[1], p[2] + q[2])


def inside(x, k, l, m, p, q):
    a, b, c = x
    s = a + b - c
    return 0 <= a <= k and 0 <= b <= l and 0 <= c <= m and p <= s <= q


def on_boundary(x, k, l, m, p, q):
    a, b, c = x
    s = a + b - c
    return (
        a == 0 or a == k or b == 0 or b == l or c == 0 or c == m
        or s == p or s == q
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
        "s=p": [(Z, BC, AC), (Z, AC, AMB)],
    }


def facet_specs(k, l, m, p, q):
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
        ("s=p", lambda x: x[0] + x[1] - x[2] == p, lower["s=p"]),
        ("s=q", lambda x: x[0] + x[1] - x[2] == q, rev("s=p")),
    ]


def boundary_vertices(k, l, m, p, q):
    out = []
    for a in range(k + 1):
        for b in range(l + 1):
            for c in range(m + 1):
                x = (a, b, c)
                if inside(x, k, l, m, p, q) and on_boundary(x, k, l, m, p, q):
                    out.append(x)
    return out


def boundary_triangulation(k, l, m, p, q):
    verts = boundary_vertices(k, l, m, p, q)
    vertset = set(verts)
    tris = []
    seen = set()
    for _, is_on_facet, templates in facet_specs(k, l, m, p, q):
        for x in verts:
            if not is_on_facet(x):
                continue
            for t in templates:
                tri = tuple(add(x, u) for u in t)
                if any(v not in vertset for v in tri):
                    continue
                if not all(inside(v, k, l, m, p, q) for v in tri):
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
    return tris


def vertex_set_of_triangulation(tris):
    return sorted({v for tri in tris for v in tri})


def main():
    if len(sys.argv) != 6:
        print("usage: python boundary_tri.py k l m p q")
        raise SystemExit(2)
    k, l, m, p, q = map(int, sys.argv[1:])
    tris = boundary_triangulation(k, l, m, p, q)
    verts = vertex_set_of_triangulation(tris)
    print("vertices =", len(verts))
    print("triangles =", len(tris))
    print(tris)


if __name__ == "__main__":
    main()
