# undented

Enumerate, solve, and prove all prime neoplatonic solids.

A neoplatonic solid is an undented polyhedron with equilateral triangle
faces, meeting at most six to a vertex.  It is prime if it has no
degree-3 vertex.

## Quick start

    make all VMAX=12       # everything through v=12 (~1 minute)

For large v, clone to a compute server:

    ssh doob 'cd undented && git pull && nohup make all VMAX=50 JOBS=96 > run/logs/all.log 2>&1 &'
    ssh doob 'cd undented && make status'

Requires: C compiler, GNU parallel, python3, numpy, scipy, LAPACK.

## Pipeline

    make seeds   VMAX=N     fullerene-dual seeds via buckygen
    make primes  VMAX=N     grow by recurrence: prime(v) = grow(prime(v-1)) ∪ seed(v)
    make solve   VMAX=N     Euclidean embeddings via homotopy from ideal
    make prove   VMAX=N     existence proofs (IEEE 754 + error bounds)
    make status             show results

## Layout

    src/                    our code
      grow_step.c             vertex insertion → offspring names
      clers_name.c            canonical CLERS naming
      neoeuc_c.c              homotopy solver (ideal → Euclidean)
      newton_polish.c         Newton-refine to machine precision
      prove_float.py          existence prover (float + error bounds)
      bin2obj.py              solver binary → OBJ converter
      plantri_to_poly         planar_code → face-list converter
    third_party/buckygen/   Brinkmann's fullerene generator (GPL)
    data/                   generated data (gitignored)
    run/                    scratch and logs (gitignored)
    PROOF_METHOD.txt        why the float prover is rigorous
