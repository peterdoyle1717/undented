# undented

Enumerate, solve, and prove all prime neoplatonic solids.

A neoplatonic solid is an undented polyhedron with equilateral triangle
faces, meeting at most six to a vertex.  It is prime if it has no
degree-3 vertex.

## Quick start

    make all VMAX=12       # everything through v=12 (~1 minute)

For large v, clone to a compute server:

    ssh server 'cd undented && git pull && nohup make all VMAX=50 JOBS=96 > run/logs/all.log 2>&1 &'
    ssh server 'cd undented && make status'

Requires: C compiler, GNU parallel, python3, numpy, scipy.

## Pipeline

    make seeds   VMAX=N     fullerene-dual seeds via buckygen
    make primes  VMAX=N     grow by recurrence: prime(v) = grow(prime(v-1)) ∪ seed(v)
    make solve   VMAX=N     Euclidean embeddings via homotopy from ideal
    make prove   VMAX=N     existence proofs (IEEE 754 + error bounds)
    make status             show results

## What the prover checks

0. **Undented** — every vertex turning sum > 0
0. **Embedded** — no triangle-triangle intersections
1. **3|V| ≥ |E|**
2. **Collision distance > 0** — non-adjacent simplices separated
3. **σ_min > 0 and ρ < σ_min²/(16√E)** — rigidity + near solution
4. **Perturbation bound < CD/√V** — solution stays embedded

See PROOF_METHOD.txt for details.

## Layout

    src/                    source code
      clers.c                 CLERS tool: decode, encode, name, canonical
      clers.py                CLERS decoder (Python, for scripting)
      grow_step.c             vertex insertion → offspring names
      neoeuc_c.c              homotopy solver (ideal → Euclidean, with polish)
      prove_float.py          existence prover (float + IEEE 754 error bounds)
      plantri_to_poly         planar_code → face-list converter
    third_party/buckygen/   Brinkmann's fullerene generator (GPL)
    data/                   generated data (gitignored)
    run/                    scratch and logs (gitignored)
    PROOF_METHOD.txt        why the float prover is rigorous
