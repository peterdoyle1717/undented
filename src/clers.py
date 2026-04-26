#!/usr/bin/env python3
"""
CLERS encoding and decoding of triangulations of the sphere.

A triangulation is a list of triples (tuples of 3 ints), vertices 1..n.
Each oriented edge (a,b) belongs to exactly one face (a,b,c).

CLERS alphabet: E=leaf, A/B=single child, C=handle, D=binary split.
One character per triangle.
"""

from collections import deque


# --- edge map ---

def build_edge_map(poly):
    """Dict mapping oriented edge (a,b) -> opposite vertex c."""
    em = {}
    for a, b, c in poly:
        em[(a, b)] = c
        em[(b, c)] = a
        em[(c, a)] = b
    return em


# --- encode (iterative) ---

def encode(poly, start=None):
    em = build_edge_map(poly)
    if start is None:
        a0, b0 = poly[0][0], poly[0][1]
    else:
        a0, b0 = start

    vset = set()
    tree = set()
    code = []

    vset.add(a0)
    vset.add(b0)
    tree.add((b0, a0))

    # explicit stack replaces recursion; process edges in LIFO order
    stack = [(a0, b0)]

    while stack:
        x, y = stack.pop()
        z = em[(x, y)]

        if z not in vset:
            code.append("C")
            tree.add((z, y))
            vset.add(z)
            stack.append((x, z))          # tail call

        elif (y, z) in tree and (z, x) in tree:
            code.append("E")

        elif (z, x) in tree:
            code.append("A")
            stack.append((z, y))          # tail call

        elif (y, z) in tree:
            code.append("B")
            stack.append((x, z))          # tail call

        else:
            code.append("D")
            # push second call first (LIFO), so first call runs next
            stack.append((z, y))
            stack.append((x, z))

    return "".join(code)


# --- decode (iterative) ---

# Stack frame states for the iterative decoder:
# Each frame = (tile, a, b, c, phase)
# phase tracks which recursive call we're waiting on.
#   A,B,C: phase 0 = push child, phase 1 = fixup + emit
#   D:     phase 0 = push first child, phase 1 = push second, phase 2 = merge + emit
#   E:     handled immediately (no children)

_VALID = frozenset("EABCD")


def decode(recipe, verify=False):
    n = len(recipe)
    for i, c in enumerate(recipe):
        if c not in _VALID:
            raise ValueError(
                f"invalid CLERS character {c!r} at position {i} "
                f"(valid alphabet: {''.join(sorted(_VALID))})"
            )
    parent = {}       # union-find, inline for speed
    triangles = []    # appended in order, reversed at end

    # inline union-find
    def find(x):
        while x in parent:
            p = parent[x]
            if p == x:
                break
            # path splitting: point x to grandparent
            gp = parent.get(p, p)
            parent[x] = gp
            x = gp
        return x

    def link(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    ptr = 0
    stack = []        # (tile, a, b, c, phase)
    dq_stack = []     # stack of deque results

    # seed
    V = 2
    v0, v1 = 1, 2

    # read first tile and push
    tile = recipe[ptr]; ptr += 1
    V += 1; c0 = V
    stack.append((tile, v0, v1, c0, 0))

    while stack:
        tile, a, b, c, phase = stack[-1]

        if tile == "E":
            stack.pop()
            dq_stack.append(deque([a, c, b]))
            triangles.append((a, b, c))

        elif tile == "A":
            if phase == 0:
                # push child chew(c, b)
                stack[-1] = (tile, a, b, c, 1)
                t = recipe[ptr]; ptr += 1
                V += 1
                stack.append((t, c, b, V, 0))
            else:
                stack.pop()
                dq = dq_stack[-1]      # result from child, reuse in place
                dq.appendleft(a)
                triangles.append((a, b, c))

        elif tile == "B":
            if phase == 0:
                stack[-1] = (tile, a, b, c, 1)
                t = recipe[ptr]; ptr += 1
                V += 1
                stack.append((t, a, c, V, 0))
            else:
                stack.pop()
                dq = dq_stack[-1]
                dq.append(b)
                triangles.append((a, b, c))

        elif tile == "C":
            if phase == 0:
                stack[-1] = (tile, a, b, c, 1)
                t = recipe[ptr]; ptr += 1
                V += 1
                stack.append((t, a, c, V, 0))
            else:
                stack.pop()
                dq = dq_stack[-1]
                dq.pop()              # drop trailing c
                d = dq.pop()          # penultimate
                dq.append(b)
                link(d, b)
                triangles.append((a, b, c))

        elif tile == "D":
            if phase == 0:
                # push first child chew(a, c)
                stack[-1] = (tile, a, b, c, 1)
                t = recipe[ptr]; ptr += 1
                V += 1
                stack.append((t, a, c, V, 0))
            elif phase == 1:
                # push second child chew(c, b)
                stack[-1] = (tile, a, b, c, 2)
                t = recipe[ptr]; ptr += 1
                V += 1
                stack.append((t, c, b, V, 0))
            else:
                stack.pop()
                dq2 = dq_stack.pop()
                dq = dq_stack[-1]
                dq.pop()              # drop shared c
                dq.extend(dq2)
                triangles.append((a, b, c))

    assert ptr == n, f"leftover input at position {ptr}"
    assert len(dq_stack) == 1

    # resolve identifications and relabel
    triangles.reverse()
    poly = [(find(a), find(b), find(c)) for a, b, c in triangles]
    poly = _regulate(poly)

    if verify:
        assert encode(poly) == recipe, "round-trip failed"

    return poly


# --- canonical naming ---

def _regulate(poly):
    """Relabel vertices 1..n in order of first appearance."""
    label = {}
    k = 0
    out = []
    for a, b, c in poly:
        if a not in label:
            k += 1; label[a] = k
        if b not in label:
            k += 1; label[b] = k
        if c not in label:
            k += 1; label[c] = k
        out.append((label[a], label[b], label[c]))
    return out


def min_degree_edges(poly):
    """Directed edges (a,b) where a has minimum degree."""
    deg = {}
    for a, b, c in poly:
        for v in (a, b, c):
            deg[v] = deg.get(v, 0) + 1
    min_deg = min(deg.values())
    em = build_edge_map(poly)
    return [e for e in em if deg[e[0]] == min_deg]


def reverse_poly(poly):
    """Reverse orientation: each triangle (a,b,c) becomes (c,b,a)."""
    return [(c, b, a) for a, b, c in poly]


def official_name(poly):
    """Lex-first CLERS string, starting only from minimum-degree vertices."""
    return min(encode(poly, start=e) for e in min_degree_edges(poly))


def official_unoriented_name(poly):
    """Lex-first CLERS string across both orientations."""
    return min(official_name(poly), official_name(reverse_poly(poly)))
