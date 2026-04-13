#!/usr/bin/env python3
"""Convert neoeuc_c binary output to OBJ files.

Usage: python3 bin2obj.py prime/N.txt euc/bin/N.bin outdir/
"""
import sys, os, struct
import numpy as np

prime_file, bin_file, outdir = sys.argv[1:4]
os.makedirs(outdir, exist_ok=True)

import gzip
opener = gzip.open if prime_file.endswith('.gz') else open
names = [l.strip() for l in opener(prime_file, 'rt') if l.strip()]

# Decode face lists
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from clers import decode

data = open(bin_file, 'rb').read()
offset = 0
count = 0
failed = 0
for name in names:
    poly = decode(name)
    nv = max(max(f) for f in poly)
    nbytes = nv * 3 * 8
    if offset + nbytes > len(data):
        break
    coords = np.frombuffer(data[offset:offset+nbytes], dtype=np.float64).reshape(nv, 3)
    offset += nbytes
    if np.isnan(coords[0, 0]):
        # Write failure marker so the directory tells the full story
        open(os.path.join(outdir, name + '.failed'), 'w').close()
        failed += 1
        continue
    objpath = os.path.join(outdir, name + '.obj')
    with open(objpath, 'w') as f:
        for v in coords:
            f.write(f'v {v[0]:.15g} {v[1]:.15g} {v[2]:.15g}\n')
        for a, b, c in poly:
            f.write(f'f {a} {b} {c}\n')
    count += 1

v_label = os.path.basename(bin_file).replace(".bin", "")
print(f'v={v_label}: {count} OBJs, {failed} failed')
