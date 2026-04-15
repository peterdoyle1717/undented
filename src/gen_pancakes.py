#!/usr/bin/env python3
"""
gen_pancakes.py — generate CLERS codes of all prime pancakes up to v vertices.

A prime pancake is the double cover of a convex region in the triangular
(Eisenstein) lattice, with all boundary angles 120° or 180°.

The region is a convex hexagon with side lengths (a,b,c,d,e,f) along the
six lattice directions (cycling 60° each), where opposite sides satisfy:
  a+b = d+e, b+c = e+f, c+d = f+a  (closure condition)
which simplifies to a+b = d+e, c = a, f = b (so a=c, b=f, d+e = a+b).
Wait, let me re-derive...

A convex hex in the triangular lattice has 6 sides with directions
rotating 60° each. If side lengths are s1..s6, closure gives:
  s1 + s2*w + s3*w^2 + s4*w^3 + s5*w^4 + s6*w^5 = 0
where w = e^(i*pi/3). This gives: s1=s4, s2=s5, s3=s6 (opposite sides equal)
Wait, that's only for regular hexagons. For general convex hex with 120° angles:

Directions: d0=1, d1=w, d2=w^2, d3=-1, d4=-w, d5=-w^2 where w=e^(i*pi/3).
Closure: sum(s_k * d_k) = 0. Since d3=-d0, d4=-d1, d5=-d2:
  (s1-s4)*d0 + (s2-s5)*d1 + (s3-s6)*d2 = 0
Since d0, d1, d2 are linearly independent over R (well, d0+d1+d2 = ... hmm)
Actually d0=1, d1=(1+i√3)/2, d2=(-1+i√3)/2. Real part: s1-s4 + (s2-s5)/2 - (s3-s6)/2 = 0.
Imag part: (s2-s5)√3/2 + (s3-s6)√3/2 = 0, so s2-s5 = -(s3-s6), i.e. s2+s6 = s3+s5.
And from real: s1-s4 + (s2-s5-s3+s6)/2 = 0, combined with s2+s6=s3+s5:
  s2-s5 = s6-s3, so s1-s4 + (s6-s3-s3+s6)/2 = s1-s4+s6-s3 = 0, so s4 = s1+s6-s3.

So the constraints are:
  s5 = s2 + s6 - s3
  s4 = s1 + s6 - s3
with s4, s5 >= 0.

Actually I think the standard parameterization is simpler. A convex hex region
in the triangular lattice is determined by three parameters (a, b, c) where:
- Walk: a steps in dir 0, b steps in dir 1, c steps in dir 2,
        then back: (b+c-a) steps in dir 3, ... no this isn't right either.

Let me just use coordinates directly. The triangular lattice has basis
e1 = (1, 0), e2 = (1/2, sqrt(3)/2). A lattice point is p = m*e1 + n*e2.

A convex region: enumerate all convex subsets of lattice points up to
size ~v/2, take double cover, compute face list, encode CLERS.

For v up to 50, the region has at most ~25 interior + boundary points.
Brute force: for each convex lattice polygon with k points where
2*(k - boundary) + boundary <= 50, generate and encode.

Actually, let me parameterize hexagons directly.
"""

import sys, os, subprocess, math
from itertools import product

def eisenstein_hex_vertices(a, b, c):
    """Generate vertices of convex hex region with side lengths a,b,c,a,b,c.

    The hexagon has vertices at:
      P0 = (0, 0)
      P1 = P0 + a * dir0
      P2 = P1 + b * dir1
      P3 = P2 + c * dir2
      P4 = P3 + a * dir3  (= P3 - a * dir0)
      P5 = P4 + b * dir4  (= P4 - b * dir1)
      P6 = P5 + c * dir5  (should = P0)

    where dir0 = (1,0), dir1 = (1/2, sqrt(3)/2), dir2 = (-1/2, sqrt(3)/2).
    Only valid if a,b,c >= 1.
    For a=0 or b=0 or c=0 we get degenerate (triangle or rhombus).
    """
    # Lattice directions (in Cartesian)
    dirs = [(1, 0), (0.5, math.sqrt(3)/2), (-0.5, math.sqrt(3)/2),
            (-1, 0), (-0.5, -math.sqrt(3)/2), (0.5, -math.sqrt(3)/2)]

    # Walk the boundary
    sides = [a, b, c, a, b, c]
    boundary = []
    x, y = 0.0, 0.0
    for side_idx in range(6):
        dx, dy = dirs[side_idx]
        for step in range(sides[side_idx]):
            boundary.append((round(2*x)/2, round(2*y*2/math.sqrt(3))/2*math.sqrt(3)/2))
            x += dx
            y += dy
    return boundary

def lattice_points_in_hex_general(s1, s2, s3, s4, s5, s6):
    """All triangular lattice points in hex with sides s1..s6.

    In axial coords (q,r), lattice directions are:
      dir0=(1,0), dir1=(0,1), dir2=(-1,1), dir3=(-1,0), dir4=(0,-1), dir5=(1,-1).

    Hex vertices:
      P0 = (0,0)
      P1 = P0 + s1*dir0 = (s1, 0)
      P2 = P1 + s2*dir1 = (s1, s2)
      P3 = P2 + s3*dir2 = (s1-s3, s2+s3)
      P4 = P3 + s4*dir3 = (s1-s3-s4, s2+s3)
      P5 = P4 + s5*dir4 = (s1-s3-s4, s2+s3-s5)
      P6 = P5 + s6*dir5 = (s1-s3-s4+s6, s2+s3-s5-s6) = (0,0) ✓
    """
    vertices = [
        (0, 0),
        (s1, 0),
        (s1, s2),
        (s1-s3, s2+s3),
        (s1-s3-s4, s2+s3),
        (s1-s3-s4, s2+s3-s5),
    ]
    # Bounding box
    qs = [v[0] for v in vertices]
    rs = [v[1] for v in vertices]
    qmin, qmax = min(qs), max(qs)
    rmin, rmax = min(rs), max(rs)

    def in_hex(q, r):
        n = len(vertices)
        for i in range(n):
            x1, y1 = vertices[i]
            x2, y2 = vertices[(i+1) % n]
            cross = (x2-x1)*(r-y1) - (y2-y1)*(q-x1)
            if cross < -1e-9:
                return False
        return True

    points = []
    for q in range(qmin, qmax+1):
        for r in range(rmin, rmax+1):
            if in_hex(q, r):
                points.append((q, r))
    return points, vertices

def lattice_points_in_hex(a, b, c):
    """All triangular lattice points inside or on the hex with sides a,b,c.

    Use axial coordinates (q, r) where the point is q*e1 + r*e2,
    e1 = (1,0), e2 = (1/2, sqrt(3)/2).
    """
    # The hex region: walk boundary, find bounding box, test containment
    # Easier: the hex with sides a,b,c,a,b,c has constraints in axial coords:
    #   0 <= q <= a+b
    #   0 <= r <= b+c
    #   q - r <= a      (from the top-right constraint)
    #   r - q <= c      (from the bottom-left constraint)
    #   q + (b+c-r) <= a+b  =>  q <= r+a  =>  q-r <= a  (same)
    # Let me think again...

    # Place P0 at origin. The hex vertices in axial (q,r):
    # P0 = (0,0)
    # P1 = (a, 0)       [a steps in dir0 = e1]
    # P2 = (a+b, b)     [b steps in dir1 = e2... wait, dir1 = e1+e2? No.]

    # e1 = (1,0) in axial = dir0
    # e2 = (0,1) in axial. In Cartesian: (1/2, sqrt(3)/2) = dir1. Good.
    # dir2 in Cartesian is (-1/2, sqrt(3)/2) = e2 - e1 in axial = (-1, 1).

    # So in axial coords:
    # dir0 = (1,0), dir1 = (0,1), dir2 = (-1,1),
    # dir3 = (-1,0), dir4 = (0,-1), dir5 = (1,-1).

    # Hex vertices in axial:
    # P0 = (0, 0)
    # P1 = (a, 0)
    # P2 = (a, b)       [b steps in dir1=(0,1)]
    # P3 = (a-c, b+c)   [c steps in dir2=(-1,1)]
    # P4 = (-c, b+c)    [a steps in dir3=(-1,0)]
    # P5 = (-c, c)      [b steps in dir4=(0,-1)]
    # Back to P0: c steps in dir5=(1,-1): (-c+c, c-c) = (0,0). ✓

    # The convex hull of these 6 points. Interior lattice points satisfy:
    # Half-plane constraints from each edge.
    # But simpler: in axial coords, the constraints are:
    #   r >= 0           (below P0-P1)
    #   q <= a           (right of P1-P2... no, P1=(a,0), P2=(a,b), so q<=a? Yes for this edge)
    # Actually let me just enumerate the bounding box and check point-in-polygon.

    vertices = [(0,0), (a,0), (a,b), (a-c,b+c), (-c,b+c), (-c,c)]
    # Bounding box in axial
    qs = [v[0] for v in vertices]
    rs = [v[1] for v in vertices]
    qmin, qmax = min(qs), max(qs)
    rmin, rmax = min(rs), max(rs)

    # Point-in-convex-polygon test using cross products
    def in_hex(q, r):
        # Test if (q,r) is inside or on the convex hull of vertices
        n = len(vertices)
        for i in range(n):
            x1, y1 = vertices[i]
            x2, y2 = vertices[(i+1) % n]
            cross = (x2-x1)*(r-y1) - (y2-y1)*(q-x1)
            if cross < -1e-9:
                return False
        return True

    points = []
    for q in range(qmin, qmax+1):
        for r in range(rmin, rmax+1):
            if in_hex(q, r):
                points.append((q, r))
    return points, vertices

def double_cover_faces(points, vertices_hex):
    """Build face list of double cover of a triangulated hex region.

    The region is triangulated by the triangular lattice.
    Interior points are doubled (top/bottom copy).
    Boundary points are shared.
    """
    pts = set(points)
    # Identify boundary points
    boundary = set()
    # A point is on the boundary if it's missing at least one of its 6 neighbors
    nbr_dirs = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
    for p in points:
        for d in nbr_dirs:
            if (p[0]+d[0], p[1]+d[1]) not in pts:
                boundary.add(p)
                break
    interior = pts - boundary

    # Vertex numbering: boundary points get one index,
    # interior points get two (top=original, bottom=new)
    idx = {}
    n = 0
    for p in sorted(points):
        n += 1
        idx[('T', p)] = n  # top copy
        if p in interior:
            n += 1
            idx[('B', p)] = n  # bottom copy
        else:
            idx[('B', p)] = idx[('T', p)]  # boundary: shared

    nv = n
    # Triangles: for each upward triangle (q,r)-(q+1,r)-(q,r+1)
    # and each downward triangle (q+1,r)-(q+1,r+1)-(q,r+1)
    faces_top = []
    faces_bot = []
    qs = [p[0] for p in points]
    rs = [p[1] for p in points]
    qmin, qmax = min(qs), max(qs)
    rmin, rmax = min(rs), max(rs)
    for q in range(qmin-1, qmax+1):
        for r in range(rmin-1, rmax+1):
            # Upward: (q,r), (q+1,r), (q,r+1) — all three must be in pts
            if (q,r) in pts and (q+1,r) in pts and (q,r+1) in pts:
                faces_top.append((idx[('T',(q,r))], idx[('T',(q+1,r))], idx[('T',(q,r+1))]))
                faces_bot.append((idx[('B',(q,r))], idx[('B',(q,r+1))], idx[('B',(q+1,r))]))
            # Downward: (q+1,r), (q+1,r+1), (q,r+1) — all three must be in pts
            if (q+1,r) in pts and (q+1,r+1) in pts and (q,r+1) in pts:
                faces_top.append((idx[('T',(q+1,r))], idx[('T',(q+1,r+1))], idx[('T',(q,r+1))]))
                faces_bot.append((idx[('B',(q+1,r))], idx[('B',(q,r+1))], idx[('B',(q+1,r+1))]))

    faces = faces_top + faces_bot
    return nv, faces

def faces_to_clers(nv, faces, clers_bin):
    """Convert face list to canonical CLERS code using bin/clers."""
    facelist = ';'.join(f'{a},{b},{c}' for a,b,c in faces)
    # Use clers encode then clers canonical
    result = subprocess.run([clers_bin, 'name'],
                          input=facelist, capture_output=True, text=True)
    return result.stdout.strip()

def main():
    import argparse
    ap = argparse.ArgumentParser(description='Generate prime pancake CLERS codes')
    ap.add_argument('--vmax', type=int, default=50)
    ap.add_argument('--clers', default='bin/clers', help='path to clers binary')
    args = ap.parse_args()

    pancakes = set()
    # General convex hex in triangular lattice: sides (s1,s2,s3,s4,s5,s6)
    # with closure: s4 = s1+s6-s3, s5 = s2+s6-s3.
    # All sides >= 1 (0-length side creates degree-2 vertex).
    # Free parameters: s1, s2, s3, s6.

    for s1 in range(1, args.vmax):
        for s2 in range(1, args.vmax):
            for s3 in range(1, args.vmax):
                for s6 in range(1, args.vmax):
                    s4 = s1 + s6 - s3
                    s5 = s2 + s6 - s3
                    if s4 < 1 or s5 < 1:
                        continue

                    # Quick size estimate: area ~ s1*s2 + s2*s3 + ...
                    # More precisely, total lattice points ~ area + perimeter/2 + 1
                    perim = s1+s2+s3+s4+s5+s6
                    if perim > 2*args.vmax:
                        break

                    try:
                        points, hex_verts = lattice_points_in_hex_general(s1,s2,s3,s4,s5,s6)
                    except:
                        continue

                    npts = len(points)
                    if npts < 4:
                        continue

                    pts_set = set(points)
                    nbr_dirs = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
                    nboundary = sum(1 for p in points
                                    if any((p[0]+d[0],p[1]+d[1]) not in pts_set for d in nbr_dirs))
                    ninterior = npts - nboundary
                    nv = 2 * ninterior + nboundary

                    if nv < 4 or nv > args.vmax:
                        continue

                    try:
                        nv_dc, faces = double_cover_faces(points, hex_verts)
                    except:
                        continue

                    if nv_dc > args.vmax or nv_dc < 4:
                        continue
                    if len(faces) != 2*nv_dc - 4:
                        continue  # not a valid sphere

                    code = faces_to_clers(nv_dc, faces, args.clers)
                    if code:
                        pancakes.add(code)

    for code in sorted(pancakes, key=lambda c: (len(c), c)):
        v = (len(code) + 4) // 2
        print(f'{code}')

if __name__ == '__main__':
    main()
