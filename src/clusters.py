"""Identify face clusters: connected components of coplanar adjacent triangles.

For each flop OBJ:
  - Compute per-triangle normal
  - Build adjacency: triangles sharing an edge AND coplanar (parallel normals, same plane)
  - Connected components = face clusters
  - For each cluster: report size (n triangles), interior degree-6 vertex count
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
    nf = len(faces)
    nv = len(verts)
    # Per-face normal and plane offset
    norms = np.zeros((nf, 3))
    offs = np.zeros(nf)
    for i, f in enumerate(faces):
        v1 = verts[f[1]] - verts[f[0]]
        v2 = verts[f[2]] - verts[f[0]]
        n = np.cross(v1, v2)
        nm = np.linalg.norm(n)
        if nm > 1e-12:
            n /= nm
        norms[i] = n
        offs[i] = -np.dot(n, verts[f[0]])

    # Edge → faces sharing it
    edge_to_faces = {}
    for i, f in enumerate(faces):
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            e = (min(u,v), max(u,v))
            edge_to_faces.setdefault(e, []).append(i)

    # Union-find for clusters
    parent = list(range(nf))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb

    for e, fs in edge_to_faces.items():
        if len(fs) != 2: continue
        a, b = fs
        # Coplanar test: parallel normals AND same offset
        # |n_a × n_b| ≈ 0 AND |offs_a - offs_b * sign(n_a·n_b)| ≈ 0
        # Equivalently: the third vertex of one face lies on the plane of the other
        # Get the vertex of face b that is NOT on the shared edge
        opp_b = [v for v in faces[b] if v not in (e[0], e[1])][0]
        d = abs(norms[a] @ verts[opp_b] + offs[a])
        if d < eps:
            union(a, b)

    # Build clusters
    clusters = {}
    for i in range(nf):
        r = find(i)
        clusters.setdefault(r, []).append(i)

    return list(clusters.values())

def cluster_info(verts, faces, cluster):
    """For a cluster of face indices, count vertices and interior deg-6."""
    cluster_faces = [faces[i] for i in cluster]
    cluster_verts = set()
    for f in cluster_faces:
        cluster_verts.update(f)

    # Vertex degree within cluster
    deg = {v: 0 for v in cluster_verts}
    nbrs = {v: set() for v in cluster_verts}
    for f in cluster_faces:
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            nbrs[u].add(v); nbrs[v].add(u)
    for v in cluster_verts:
        deg[v] = len(nbrs[v])

    # Boundary of cluster: edges that appear in only one cluster face
    edge_count = {}
    for f in cluster_faces:
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            e = (min(u,v), max(u,v))
            edge_count[e] = edge_count.get(e, 0) + 1
    boundary_verts = set()
    for e, c in edge_count.items():
        if c == 1:
            boundary_verts.add(e[0]); boundary_verts.add(e[1])

    interior_verts = cluster_verts - boundary_verts
    interior_deg6 = sum(1 for v in interior_verts if deg[v] == 6)

    return len(cluster_faces), len(cluster_verts), len(boundary_verts), len(interior_verts), interior_deg6

if __name__ == '__main__':
    for path in sys.argv[1:]:
        verts, faces = read_obj(path)
        clusters = face_clusters(verts, faces)
        clusters.sort(key=len, reverse=True)
        infos = [cluster_info(verts, faces, c) for c in clusters]
        # Show clusters with > 1 face
        nontriv = [(len(c), info) for c, info in zip(clusters, infos) if len(c) > 1]
        name = os.path.basename(path).replace('.obj','')
        v = (len(name) + 4) // 2
        if nontriv:
            print(f'v={v} {name}:')
            for nf_c, info in nontriv:
                nf2, nv2, nb, ni, ni6 = info
                print(f'  cluster: {nf2} tris, {nv2} verts ({nb} bdry, {ni} int, {ni6} int-deg6)')
        else:
            print(f'v={v} {name}: no nontrivial clusters')
