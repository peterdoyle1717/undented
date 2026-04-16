"""Build cluster database from tame flops, look for gluings.

For each tame flop:
  - Compute clusters
  - For each cluster: signature, corner degrees, vertex count, face count

Two clusters can glue if their signatures are mirror images (or we match
the canonical form which already accounts for reflection). Corner
degrees at matched vertices must produce final degree in {4,5,6}.
"""
import sys, os
import numpy as np
sys.path.insert(0, '/Users/doyle/Library/CloudStorage/Dropbox/neo/undented/undented/src')
from cluster_sig import face_clusters, cluster_boundary_signature, read_obj

def cluster_details(verts, faces, cluster):
    """Return (sig, corner_degs_tuple, n_tris, n_verts).

    corner_degs_tuple: for each boundary vertex in cyclic order,
    its (degree_in_flop, k_in_cluster).
    """
    cluster_faces_set = set(cluster)
    cluster_verts = set()
    for fi in cluster:
        cluster_verts.update(faces[fi])

    # Full polyhedron degrees
    full_nbrs = {i: set() for i in range(len(verts))}
    for f in faces:
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            full_nbrs[u].add(v); full_nbrs[v].add(u)

    # Boundary edges of cluster: appear in exactly one cluster face
    edge_count = {}
    for fi in cluster:
        f = faces[fi]
        for u, v in [(f[0],f[1]),(f[1],f[2]),(f[0],f[2])]:
            e = (min(u,v), max(u,v))
            edge_count[e] = edge_count.get(e, 0) + 1
    boundary_edges = [e for e, c in edge_count.items() if c == 1]
    if not boundary_edges:
        return None  # no boundary (e.g. pancake)

    # k_in_cluster at each cluster vertex = number of cluster triangles incident
    k_in_cluster = {v: 0 for v in cluster_verts}
    for fi in cluster:
        for v in faces[fi]:
            k_in_cluster[v] += 1

    # Walk boundary
    bdry_oriented = []
    for fi in cluster:
        f = faces[fi]
        for i in range(3):
            u, v = f[i], f[(i+1)%3]
            if (min(u,v), max(u,v)) in set(boundary_edges):
                bdry_oriented.append((u, v))
    nxt = {u: v for u, v in bdry_oriented}
    if len(nxt) != len(bdry_oriented):
        return None
    start = bdry_oriented[0][0]
    cycle = [start]
    cur = nxt[start]
    while cur != start:
        cycle.append(cur)
        if cur not in nxt: return None
        cur = nxt[cur]
    if len(cycle) != len(bdry_oriented):
        return None

    sig = cluster_boundary_signature(verts, faces, cluster)
    corner_info = tuple((len(full_nbrs[v]), k_in_cluster[v]) for v in cycle)
    return (sig, corner_info, len(cluster), len(cluster_verts))

def analyze_flop(name):
    """Read OBJ for a flop, return list of (sig, corner_info, ntris, nverts) per cluster."""
    import subprocess
    v = (len(name) + 4) // 2
    obj = subprocess.run(['ssh','-n','doob.dartmouth.edu',f'cat ~/undented/data/obj/{v}/{name}.obj'],
                       capture_output=True, text=True, timeout=30)
    if obj.returncode != 0:
        return None, None, None
    lines = obj.stdout.split('\n')
    verts, faces = [], []
    for line in lines:
        p = line.split()
        if not p: continue
        if p[0] == 'v':
            verts.append([float(x) for x in p[1:4]])
        elif p[0] == 'f':
            faces.append([int(x.split('/')[0]) - 1 for x in p[1:4]])
    verts = np.array(verts)
    clusters = face_clusters(verts, faces)
    clusters = [c for c in clusters if len(c) > 1]
    details = []
    for c in clusters:
        d = cluster_details(verts, faces, c)
        if d:
            details.append(d)
    return details, len(verts), len(faces)

if __name__ == '__main__':
    for name in sys.argv[1:]:
        details, nv, nf = analyze_flop(name)
        v = (len(name) + 4) // 2
        if not details:
            print(f'{name} (v={v}): no clusters')
            continue
        print(f'{name} (v={v}, F={nf}):')
        for sig, corner_info, ntris, nverts in details:
            corners = [(d, k) for d, k in corner_info if k <= 2]  # only true corners
            print(f'  sig={sig}')
            print(f'    {ntris}t/{nverts}v corners (d,k)={corner_info}')
