#!/usr/bin/env python3
"""gen_preapproved.py — emit per-v CLERS lists of all known winner floppers.

Inductively builds the set of "preapproved" (= solved-by-construction)
floppers up to --vmax:

  atomic
    pancake          gen_pancakes       (double cover of lattice polygon)
    octahedron       boundary_tri       (cannonball region; require_no_deg3)
    thick pancake    gen_thick_pancakes (planar polygon × [0,1] hip-roof)
    misc             v=14 hex antiprism (CCCCACCACACACACACAACAAAE)

  compound  (one round of cfglue, donor pool = atomic + deg-3-allowed octs)

Output: data/preapproved/<v>.txt, one CLERS per line, sorted.

Usage: python3 gen_preapproved.py [--vmax 50] [--outdir data/preapproved]
"""
from __future__ import annotations
import sys, os, argparse
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from clers import official_unoriented_name
from gen_pancakes import enumerate_pancakes, pancake_facelist, edge_lengths, canonical_sextuple
from gen_thick_pancakes import enumerate_polygons, build_thick_pancake
from gen_flop_gallery import (compute_pancakes, compute_octahedra,
                              octahedron_donors, compute_compounds_from_donors)
from cfglue import extract_bigfaces

MISC = {'CCCCACCACACACACACAACAAAE'}  # v=14 hex antiprism


def vof(clers):
    return (len(clers) + 4) // 2


def gen_thick_pancake_codes(vmax):
    seen = set()
    for k, l, r, q in enumerate_polygons(vmax):
        res = build_thick_pancake(k, l, r, q)
        if res is None:
            continue
        V, faces = res
        if V > vmax:
            continue
        if len(faces) != 2 * V - 4:
            continue
        seen.add(official_unoriented_name(faces))
    return seen


def thick_pancake_donors(vmax):
    """Yield (clers, verts, faces) for each thick pancake. Verts are
    needed by extract_bigfaces in the compound step."""
    import numpy as np
    out = []
    seen = set()
    for k, l, r, q in enumerate_polygons(vmax):
        res = build_thick_pancake(k, l, r, q)
        if res is None:
            continue
        V, faces = res
        if V > vmax or len(faces) != 2 * V - 4:
            continue
        code = official_unoriented_name(faces)
        if code in seen:
            continue
        seen.add(code)
        # Use synthetic verts since we just need a 3D embedding; thick
        # pancake builder doesn't return coords. Skip for now; donor list
        # uses octs only by default.
    return out  # empty until thick-pancake coords are wired


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vmax', type=int, default=50)
    ap.add_argument('--outdir', default=os.path.join(
        os.path.dirname(__file__), '..', 'data', 'preapproved'))
    ap.add_argument('--no-compound', action='store_true',
                    help='skip compound generation (much faster)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f'computing pancakes vmax={args.vmax} ...', file=sys.stderr)
    pancakes = compute_pancakes(args.vmax)
    print(f'  {len(pancakes)} pancake CLERSes', file=sys.stderr)

    print(f'computing octahedra (no-deg3) vmax={args.vmax} ...', file=sys.stderr)
    octas = compute_octahedra(args.vmax)
    print(f'  {len(octas)} octahedron CLERSes', file=sys.stderr)

    print(f'computing thick pancakes vmax={args.vmax} ...', file=sys.stderr)
    thicks = gen_thick_pancake_codes(args.vmax)
    print(f'  {len(thicks)} thick pancake CLERSes', file=sys.stderr)

    print(f'misc: {len(MISC)} (v=14 hex antiprism)', file=sys.stderr)

    compounds = set()
    if not args.no_compound:
        print(f'computing compounds (donor pool = octahedra deg-3-allowed) ...', file=sys.stderr)
        donors = octahedron_donors(args.vmax)
        print(f'  donor pool: {len(donors)}', file=sys.stderr)
        compounds = compute_compounds_from_donors(donors)
        # filter compounds to v <= vmax
        compounds = {c for c in compounds if vof(c) <= args.vmax}
        print(f'  {len(compounds)} compound CLERSes', file=sys.stderr)

    # Combine, dedupe, group by v
    all_codes = pancakes | octas | thicks | MISC | compounds
    by_v = defaultdict(set)
    for c in all_codes:
        by_v[vof(c)].add(c)

    # Per-class breakdown counters
    print(f'\nper-v counts (total / pancake / octahedron / thick / misc / compound):', file=sys.stderr)
    total = 0
    for v in sorted(by_v):
        codes = sorted(by_v[v])
        outpath = os.path.join(args.outdir, f'{v}.txt')
        with open(outpath, 'w') as f:
            for c in codes:
                f.write(c + '\n')
        np = sum(1 for c in codes if c in pancakes)
        no = sum(1 for c in codes if c in octas)
        nt = sum(1 for c in codes if c in thicks)
        nm = sum(1 for c in codes if c in MISC)
        nc = sum(1 for c in codes if c in compounds)
        print(f'  v={v:3d}  total={len(codes):6d}  '
              f'p={np:5d}  o={no:5d}  t={nt:5d}  m={nm}  c={nc:5d}',
              file=sys.stderr)
        total += len(codes)
    print(f'\ntotal CLERSes: {total}', file=sys.stderr)
    print(f'wrote per-v files to {args.outdir}', file=sys.stderr)


if __name__ == '__main__':
    main()
