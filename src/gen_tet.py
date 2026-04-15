#!/usr/bin/env python3
"""
gen_tet.py — generate CLERS codes of truncated tet subdivisions.

An n-subdivision of the regular tetrahedron tiles each face with n²
equilateral triangles. Truncating the 4 corners by depths t1..t4
removes the degree-3 vertices, giving a prime neoplatonic solid.

Parameters: n (subdivision level), t1,t2,t3,t4 (truncation depths).
Constraints: ti >= 1, ti+tj <= n-1 for each tet edge (i,j).

Usage: python3 gen_tet.py --vmax 50 --clers bin/clers
"""

import sys, os, subprocess, argparse
from collections import defaultdict

def tet_subdivided_faces(n):
    """Generate the triangulated surface of an n-subdivided tetrahedron.

    Vertices are parameterized by barycentric coords (a,b,c,d) with
    a+b+c+d = n, all >= 0, and at least one coord = 0 (surface).

    Each face of the tet has vertices where one barycentric coord is 0.
    """
    # Enumerate surface vertices: (a,b,c,d) with a+b+c+d=n, min coord=0
    verts = {}
    idx = 0
    for a in range(n+1):
        for b in range(n+1-a):
            for c in range(n+1-a-b):
                d = n - a - b - c
                if min(a,b,c,d) == 0:  # surface
                    idx += 1
                    verts[(a,b,c,d)] = idx

    # Generate faces: for each tet face, triangulate the n-subdivision
    # The 4 tet faces have one coord = 0:
    # Face 0: d=0 (vertices with d=0)
    # Face 1: c=0
    # Face 2: b=0
    # Face 3: a=0
    faces = []

    # Face d=0: coords (a,b,c) with a+b+c=n
    for a in range(n):
        for b in range(n-a):
            c = n - a - b
            # Upward triangle: (a,b,c,0), (a+1,b,c-1,0), (a,b+1,c-1,0)
            if c >= 1:
                v0 = verts[(a,b,c,0)]
                v1 = verts[(a+1,b,c-1,0)]
                v2 = verts[(a,b+1,c-1,0)]
                faces.append((v0, v1, v2))
            # Downward triangle: (a+1,b,c-1,0), (a+1,b+1,c-2,0), (a,b+1,c-1,0)
            if c >= 2:
                v0 = verts[(a+1,b,c-1,0)]
                v1 = verts[(a+1,b+1,c-2,0)]
                v2 = verts[(a,b+1,c-1,0)]
                faces.append((v0, v1, v2))

    # Face c=0: coords (a,b,0,d) with a+b+d=n
    for a in range(n):
        for b in range(n-a):
            d = n - a - b
            if d >= 1:
                faces.append((verts[(a,b,0,d)], verts[(a,b+1,0,d-1)], verts[(a+1,b,0,d-1)]))
            if d >= 2:
                faces.append((verts[(a,b+1,0,d-1)], verts[(a+1,b+1,0,d-2)], verts[(a+1,b,0,d-1)]))

    # Face b=0: coords (a,0,c,d) with a+c+d=n
    for a in range(n):
        for c in range(n-a):
            d = n - a - c
            if d >= 1:
                faces.append((verts[(a,0,c,d)], verts[(a+1,0,c,d-1)], verts[(a,0,c+1,d-1)]))
            if d >= 2:
                faces.append((verts[(a+1,0,c,d-1)], verts[(a+1,0,c+1,d-2)], verts[(a,0,c+1,d-1)]))

    # Face a=0: coords (0,b,c,d) with b+c+d=n
    for b in range(n):
        for c in range(n-b):
            d = n - b - c
            if d >= 1:
                faces.append((verts[(0,b,c,d)], verts[(0,b,c+1,d-1)], verts[(0,b+1,c,d-1)]))
            if d >= 2:
                faces.append((verts[(0,b,c+1,d-1)], verts[(0,b+1,c+1,d-2)], verts[(0,b+1,c,d-1)]))

    return verts, faces

def truncate(verts, faces, n, t):
    """Remove vertices within distance t[i] of tet corner i.

    Tet corners: (n,0,0,0), (0,n,0,0), (0,0,n,0), (0,0,0,n).
    Distance of (a,b,c,d) to corner i is n - coord_i.
    Vertex is within distance t[i] of corner i if coord_i > n - t[i].
    So remove vertex if ANY coord > n - t[i] for its corner.
    Actually: vertex (a,b,c,d) is near corner 0 if a > n-t[0],
    near corner 1 if b > n-t[1], etc.
    Remove if a > n-t[0] OR b > n-t[1] OR c > n-t[2] OR d > n-t[3].
    Keep if a <= n-t[0] AND b <= n-t[1] AND c <= n-t[2] AND d <= n-t[3].
    """
    # First, add cap interior vertices that aren't on the original surface
    all_verts = dict(verts)  # copy
    next_idx = max(verts.values()) + 1
    for corner in range(4):
        ti = t[corner]
        # Cap vertices: coord[corner] = n-ti
        # Need all (a,b,c,d) with a+b+c+d=n, coord[corner]=n-ti
        coords_list = []
        if corner == 0:
            for b in range(ti+1):
                for c in range(ti+1-b):
                    d = ti-b-c
                    coords_list.append((n-ti,b,c,d))
        elif corner == 1:
            for a in range(ti+1):
                for c in range(ti+1-a):
                    d = ti-a-c
                    coords_list.append((a,n-ti,c,d))
        elif corner == 2:
            for a in range(ti+1):
                for b in range(ti+1-a):
                    d = ti-a-b
                    coords_list.append((a,b,n-ti,d))
        else:
            for a in range(ti+1):
                for b in range(ti+1-a):
                    c = ti-a-b
                    coords_list.append((a,b,c,n-ti))
        for coords in coords_list:
            if coords not in all_verts:
                all_verts[coords] = next_idx
                next_idx += 1

    keep = {}
    bary = {}
    new_idx = 0
    for (a,b,c,d), old_idx in all_verts.items():
        if a <= n-t[0] and b <= n-t[1] and c <= n-t[2] and d <= n-t[3]:
            new_idx += 1
            keep[old_idx] = new_idx
            bary[new_idx] = (a,b,c,d)

    # Reverse map: barycentric -> new_idx
    bary_to_idx = {v: k for k, v in bary.items()}

    new_faces = []
    for f in faces:
        if all(v in keep for v in f):
            new_faces.append(tuple(keep[v] for v in f))

    # Add truncation cap faces
    # Corner 0 at (n,0,0,0): cap is the triangle where a = n-t[0]
    # Cap vertices have a = n-t[0], and sit on the original tet faces
    # (so one of b,c,d = 0). The cap is a triangle with coords
    # (n-t[0], b, c, d) where b+c+d = t[0], min(b,c,d) >= 0.
    for corner in range(4):
        ti = t[corner]
        # Cap vertices: coord[corner] = n-ti, sum of others = ti, on surface
        # These form a triangulated triangle of side ti
        cap_verts = {}
        if corner == 0:
            for b in range(ti+1):
                for c in range(ti+1-b):
                    d = ti - b - c
                    key = (n-ti, b, c, d)
                    if key in bary_to_idx:
                        cap_verts[(b, c)] = bary_to_idx[key]
        elif corner == 1:
            for a in range(ti+1):
                for c in range(ti+1-a):
                    d = ti - a - c
                    key = (a, n-ti, c, d)
                    if key in bary_to_idx:
                        cap_verts[(a, c)] = bary_to_idx[key]
        elif corner == 2:
            for a in range(ti+1):
                for b in range(ti+1-a):
                    d = ti - a - b
                    key = (a, b, n-ti, d)
                    if key in bary_to_idx:
                        cap_verts[(a, b)] = bary_to_idx[key]
        else:  # corner == 3
            for a in range(ti+1):
                for b in range(ti+1-a):
                    c = ti - a - b
                    key = (a, b, c, n-ti)
                    if key in bary_to_idx:
                        cap_verts[(a, b)] = bary_to_idx[key]

        # Triangulate the cap — determine winding from adjacent surface face
        # Find an edge on the cap boundary that's shared with a surface face,
        # use that to determine orientation.
        # For simplicity: find two adjacent cap vertices that share a surface face,
        # check which winding is consistent.

        # First just collect cap faces in both windings
        cap_faces_a = []
        cap_faces_b = []
        for p in range(ti):
            for q in range(ti-p):
                if (p,q) in cap_verts and (p+1,q) in cap_verts and (p,q+1) in cap_verts:
                    cap_faces_a.append((cap_verts[(p,q)], cap_verts[(p+1,q)], cap_verts[(p,q+1)]))
                    cap_faces_b.append((cap_verts[(p,q)], cap_verts[(p,q+1)], cap_verts[(p+1,q)]))
                if (p+1,q) in cap_verts and (p+1,q+1) in cap_verts and (p,q+1) in cap_verts:
                    cap_faces_a.append((cap_verts[(p+1,q)], cap_verts[(p+1,q+1)], cap_verts[(p,q+1)]))
                    cap_faces_b.append((cap_verts[(p+1,q)], cap_verts[(p,q+1)], cap_verts[(p+1,q+1)]))

        # Check which winding is consistent: for each directed edge (u,v)
        # in the surface faces, the cap should have (v,u) not (u,v).
        surface_edges = set()
        for f in new_faces:
            surface_edges.add((f[0],f[1]))
            surface_edges.add((f[1],f[2]))
            surface_edges.add((f[2],f[0]))
        cap_vset = set(cap_verts.values())
        # Check winding a
        ok_a = True
        for f in cap_faces_a:
            for i in range(3):
                e = (f[i], f[(i+1)%3])
                if e[0] in cap_vset and e[1] in cap_vset:
                    if e in surface_edges:
                        ok_a = False
                        break
            if not ok_a:
                break
        new_faces.extend(cap_faces_a if ok_a else cap_faces_b)

    return len(keep), new_faces

def faces_to_clers(nv, faces, clers_bin):
    facelist = ';'.join(f'{a},{b},{c}' for a,b,c in faces)
    try:
        result = subprocess.run([clers_bin, 'name'],
                              input=facelist, capture_output=True, text=True,
                              timeout=5)
        return result.stdout.strip() if result.returncode == 0 else ''
    except:
        return ''

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vmax', type=int, default=50)
    ap.add_argument('--clers', default='bin/clers')
    args = ap.parse_args()

    codes = set()
    for n in range(3, args.vmax):  # n >= 3 for interior vertices
        verts, faces = tet_subdivided_faces(n)
        nv_full = len(verts)
        if nv_full > 4 * args.vmax:  # way too big
            break

        for t0 in range(1, n):
            for t1 in range(1, n):
                if t0 + t1 > n - 1:
                    continue
                for t2 in range(1, n):
                    if t0 + t2 > n - 1 or t1 + t2 > n - 1:
                        continue
                    for t3 in range(1, n):
                        if t0 + t3 > n - 1 or t1 + t3 > n - 1 or t2 + t3 > n - 1:
                            continue
                        t = [t0, t1, t2, t3]
                        nv, new_faces = truncate(verts, faces, n, t)
                        if nv < 4 or nv > args.vmax:
                            continue
                        nf = len(new_faces)
                        if nf != 2 * nv - 4:
                            continue  # not a valid sphere
                        code = faces_to_clers(nv, new_faces, args.clers)
                        if code:
                            codes.add(code)

    for code in sorted(codes, key=lambda c: (len(c), c)):
        print(code)

if __name__ == '__main__':
    main()
