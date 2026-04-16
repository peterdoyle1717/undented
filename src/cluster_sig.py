"""Boundary signature of a face cluster.

For a flat triangulated region, the boundary is a polygon. Walking the
boundary, each edge is a unit vector in one of 6 lattice directions.
The signature is the cyclic sequence of these directions, normalized
to a canonical form (smallest rotation, possibly reflected).

This signature is the same for two clusters iff they're congruent in
the lattice (same shape up to rotation/reflection).
"""
import sys, os, glob
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

def face_clusters(verts, faces, eps=1e-5):
    """Return list of clusters; each cluster is a list of face indices."""
    nf = len(faces)
    norms = np.zeros((nf, 3))
    offs = np.zeros(nf)
    for i, f in enumerate(faces):
        v1 = verts[f[1]] - verts[f[0]]
        v2 = verts[f[2]] - verts[f[0]]
        n = np.cross(v1, v2)
        nm = np.linalg.norm(n)
        if nm > 1e-12: n /= nm
        norms[i] = n
        offs[i] = -np.dot(n, verts[f[0]])

    edge_to_faces = {}
    for i, f in enumerate(faces):
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            e = (min(u,v), max(u,v))
            edge_to_faces.setdefault(e, []).append(i)

    parent = list(range(nf))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb

    for e, fs in edge_to_faces.items():
        if len(fs) != 2: continue
        a, b = fs
        opp_b = [v for v in faces[b] if v not in (e[0], e[1])][0]
        d = abs(norms[a] @ verts[opp_b] + offs[a])
        if d < eps:
            union(a, b)

    clusters = {}
    for i in range(nf):
        clusters.setdefault(find(i), []).append(i)
    return list(clusters.values())

def cluster_boundary_signature(verts, faces, cluster):
    """Compute the signature of a flat cluster's boundary.

    Returns a tuple of integers — the cyclic sequence of edge directions
    in normal form (smallest rotation/reflection).
    """
    cluster_faces = [faces[i] for i in cluster]

    # Boundary edges: appear in exactly one cluster face
    edge_count = {}
    edge_face = {}
    for fi in cluster:
        f = faces[fi]
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            e = (min(u,v), max(u,v))
            edge_count[e] = edge_count.get(e, 0) + 1
            edge_face[e] = fi

    boundary_edges = [e for e, c in edge_count.items() if c == 1]
    if not boundary_edges:
        return ()

    # Build oriented boundary by walking
    # For each face, find its boundary edges and their orientation
    bdry_oriented = []
    for fi in cluster:
        f = faces[fi]
        for i in range(3):
            u, v = f[i], f[(i+1)%3]
            if (min(u,v), max(u,v)) in set(boundary_edges):
                bdry_oriented.append((u, v))

    if not bdry_oriented: return ()

    # Walk: start at any vertex on boundary, follow oriented edges
    nxt = {u: v for u, v in bdry_oriented}
    if len(nxt) != len(bdry_oriented): return ()  # not a simple cycle

    start = bdry_oriented[0][0]
    cycle = [start]
    cur = nxt[start]
    while cur != start:
        cycle.append(cur)
        if cur not in nxt: return ()
        cur = nxt[cur]

    if len(cycle) != len(bdry_oriented): return ()  # not a single cycle

    # Convert to direction sequence: 3D unit vectors
    dirs_3d = []
    for i in range(len(cycle)):
        a = cycle[i]
        b = cycle[(i+1) % len(cycle)]
        dv = verts[b] - verts[a]
        L = np.linalg.norm(dv)
        if L < 1e-9: return ()
        dirs_3d.append(dv / L)

    # Project to 2D plane of the cluster
    # Use first face's normal
    f0 = faces[cluster[0]]
    v1 = verts[f0[1]] - verts[f0[0]]
    v2 = verts[f0[2]] - verts[f0[0]]
    n = np.cross(v1, v2); n /= np.linalg.norm(n)
    # Build 2D basis in the plane
    e1 = v1 / np.linalg.norm(v1)
    e2 = np.cross(n, e1)

    # Get 2D directions
    dirs_2d = [(d @ e1, d @ e2) for d in dirs_3d]

    # Quantize to 6 lattice directions: angle * 6 / (2π) → 0..5
    quant = []
    for x, y in dirs_2d:
        ang = np.arctan2(y, x)
        idx = round(ang * 6 / (2 * np.pi)) % 6
        quant.append(idx)

    # Canonical form: try all rotations and reflections, take min
    n = len(quant)
    candidates = []
    for rot in range(n):
        rotated = quant[rot:] + quant[:rot]
        candidates.append(tuple(rotated))
        # Reflection: reverse and negate (since direction reverses)
        reflected = [(- d) % 6 for d in rotated[::-1]]
        candidates.append(tuple(reflected))
    return min(candidates)

def cluster_summary(verts, faces, cluster):
    """Return (n_tris, n_verts, n_bdry, n_int, n_int_deg6, signature)."""
    sig = cluster_boundary_signature(verts, faces, cluster)
    cluster_verts = set()
    for fi in cluster:
        cluster_verts.update(faces[fi])
    nbrs = {v: set() for v in cluster_verts}
    for fi in cluster:
        f = faces[fi]
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            nbrs[u].add(v); nbrs[v].add(u)
    edge_count = {}
    for fi in cluster:
        f = faces[fi]
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            e = (min(u,v), max(u,v))
            edge_count[e] = edge_count.get(e, 0) + 1
    bdry = set()
    for e, c in edge_count.items():
        if c == 1:
            bdry.add(e[0]); bdry.add(e[1])
    interior = cluster_verts - bdry
    int_deg6 = sum(1 for v in interior if len(nbrs[v]) == 6)
    return (len(cluster), len(cluster_verts), len(bdry), len(interior), int_deg6, sig)

if __name__ == '__main__':
    for path in sys.argv[1:]:
        verts, faces = read_obj(path)
        clusters = face_clusters(verts, faces)
        # Filter out single-face clusters (boring)
        clusters = [c for c in clusters if len(c) > 1]
        clusters.sort(key=len, reverse=True)
        name = os.path.basename(path).replace('.obj','')
        v = (len(name) + 4) // 2
        print(f'v={v} {name}: {len(clusters)} clusters')
        for c in clusters:
            nf, nv, nb, ni, ni6, sig = cluster_summary(verts, faces, c)
            print(f'  {nf}t/{nv}v ({nb}b/{ni}i, {ni6} deg6) sig={sig}')
