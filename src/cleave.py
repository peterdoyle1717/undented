"""Detect cleavage planes: 3-point planes that split a polyhedron cleanly.

For each vertex b, for each pair of neighbors a, c of b such that (a,b,c)
is not a face, check the plane through (a,b,c):
  - some vertex strictly above, some strictly below (epsilon test)
  - no edge crosses (no edge with one endpoint strictly above, other strictly below)
"""
import sys
import numpy as np

def read_obj(path):
    verts, faces = [], []
    for line in open(path):
        p = line.split()
        if not p: continue
        if p[0] == 'v':
            verts.append([float(x) for x in p[1:4]])
        elif p[0] == 'f':
            faces.append([int(x.split('/')[0]) - 1 for x in p[1:4]])
    return np.array(verts), faces

def find_cleavages(verts, faces, eps=1e-6):
    nv = len(verts)
    # Build neighbor lists
    nbrs = {i: set() for i in range(nv)}
    face_set = set()
    for f in faces:
        a, b, c = f
        nbrs[a].add(b); nbrs[b].add(a)
        nbrs[b].add(c); nbrs[c].add(b)
        nbrs[a].add(c); nbrs[c].add(a)
        face_set.add(tuple(sorted(f)))

    cleavages = []
    seen_planes = set()
    for b in range(nv):
        nb = sorted(nbrs[b])
        for i in range(len(nb)):
            for j in range(i+1, len(nb)):
                a, c = nb[i], nb[j]
                if tuple(sorted([a,b,c])) in face_set:
                    continue
                # Plane through a,b,c
                v1 = verts[a] - verts[b]
                v2 = verts[c] - verts[b]
                n = np.cross(v1, v2)
                nm = np.linalg.norm(n)
                if nm < 1e-12: continue
                n = n / nm
                d = -np.dot(n, verts[b])
                # Signed distance of each vertex
                dist = verts @ n + d
                # Classify
                above = np.where(dist > eps)[0]
                below = np.where(dist < -eps)[0]
                on_plane = np.where(np.abs(dist) <= eps)[0]
                if len(above) == 0 or len(below) == 0:
                    continue
                # Check no edge crosses
                ok = True
                for f in faces:
                    for x, y in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
                        if dist[x] > eps and dist[y] < -eps:
                            ok = False; break
                        if dist[x] < -eps and dist[y] > eps:
                            ok = False; break
                    if not ok: break
                if not ok: continue
                # Plane signature for dedup
                sig = (tuple(sorted(on_plane)), tuple(sorted(above)))
                if sig in seen_planes: continue
                seen_planes.add(sig)
                cleavages.append((len(on_plane), len(above), len(below), sorted(on_plane.tolist())))
    return cleavages

if __name__ == '__main__':
    for path in sys.argv[1:]:
        verts, faces = read_obj(path)
        c = find_cleavages(verts, faces)
        print(f'{path}: {len(c)} cleavages, vertex counts {[(x[0],x[1],x[2]) for x in c[:3]]}')
