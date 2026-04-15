#!/usr/bin/env python3
"""
gen_pancakes.py — generate CLERS codes of all prime pancakes up to v vertices.

A pancake is the double cover of a truncated triangle in the triangular lattice.

Size-n triangle: vertices (i,j) with 0 <= i, j, i+j <= n.
Truncate corners by a, b, c >= 1 with a+b, b+c, c+a <= n-1.
  Corner at (0,0):     remove (i,j) with i+j < a
  Corner at (n,0):     remove (i,j) with i > n-b  (equiv: n-i < b)
  Corner at (0,n):     remove (i,j) with j > n-c  (equiv: n-j < c)

Double cover: two copies of each interior vertex, one copy of each
boundary vertex. Triangles from the lattice, mirrored for the second copy.
"""
import sys, subprocess, argparse

def gen_pancake(n, a, b, c, clers_bin):
    """Generate CLERS code for pancake with given parameters.
    Returns code string or '' on failure."""

    # Enumerate vertices in the truncated region
    pts = []
    for i in range(n+1):
        for j in range(n+1-i):
            if i + j < a:       continue  # corner (0,0) truncated
            if i > n - b:       continue  # corner (n,0) truncated
            if j > n - c:       continue  # corner (0,n) truncated
            pts.append((i, j))

    pts_set = set(pts)
    if len(pts) < 4:
        return ''

    # Classify: boundary vs interior
    nbr_dirs = [(1,0), (0,1), (-1,1), (-1,0), (0,-1), (1,-1)]
    boundary = set()
    for p in pts:
        for d in nbr_dirs:
            if (p[0]+d[0], p[1]+d[1]) not in pts_set:
                boundary.add(p)
                break
    interior = pts_set - boundary

    # Number vertices: boundary gets one index, interior gets two (top/bottom)
    idx = {}
    n_vert = 0
    for p in sorted(pts):
        n_vert += 1
        idx[('T', p)] = n_vert
        if p in interior:
            n_vert += 1
            idx[('B', p)] = n_vert
        else:
            idx[('B', p)] = idx[('T', p)]  # boundary: shared

    # Generate triangles
    faces = []
    imin = min(p[0] for p in pts) - 1
    imax = max(p[0] for p in pts)
    jmin = min(p[1] for p in pts) - 1
    jmax = max(p[1] for p in pts)
    for i in range(imin, imax + 1):
        for j in range(jmin, jmax + 1):
            # Upward: (i,j), (i+1,j), (i,j+1)
            if (i,j) in pts_set and (i+1,j) in pts_set and (i,j+1) in pts_set:
                faces.append((idx[('T',(i,j))], idx[('T',(i+1,j))], idx[('T',(i,j+1))]))
                faces.append((idx[('B',(i,j))], idx[('B',(i,j+1))], idx[('B',(i+1,j))]))
            # Downward: (i+1,j), (i+1,j+1), (i,j+1)
            if (i+1,j) in pts_set and (i+1,j+1) in pts_set and (i,j+1) in pts_set:
                faces.append((idx[('T',(i+1,j))], idx[('T',(i+1,j+1))], idx[('T',(i,j+1))]))
                faces.append((idx[('B',(i+1,j))], idx[('B',(i,j+1))], idx[('B',(i+1,j+1))]))

    if len(faces) != 2 * n_vert - 4:
        return ''

    facelist = ';'.join(f'{a},{b},{c}' for a,b,c in faces)
    try:
        result = subprocess.run([clers_bin, 'name'],
                              input=facelist, capture_output=True, text=True, timeout=5)
        return result.stdout.strip() if result.returncode == 0 else ''
    except:
        return ''

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vmax', type=int, default=50)
    ap.add_argument('--clers', default='bin/clers')
    args = ap.parse_args()

    codes = set()
    for n in range(3, 4*args.vmax):  # triangle side (large n with deep truncation can give small v)
        if (n+1)*(n+2)//2 - 3 > 2 * args.vmax and n > args.vmax:
            break  # even without truncation, way too big, and n itself is huge
        for a in range(1, n):
            for b in range(1, n):
                if a + b > n - 1: break
                for c in range(1, n):
                    if b + c > n - 1 or a + c > n - 1: break
                    # Count vertices quickly
                    npts = 0
                    for i in range(n+1):
                        for j in range(n+1-i):
                            if i+j < a or i > n-b or j > n-c: continue
                            npts += 1
                    # rough upper bound on v in double cover
                    if npts > args.vmax: continue
                    if npts < 4: continue

                    code = gen_pancake(n, a, b, c, args.clers)
                    if code:
                        codes.add(code)

    for code in sorted(codes, key=lambda c: (len(c), c)):
        print(code)

if __name__ == '__main__':
    main()
