# undented

Enumerate, solve, and prove all prime neoplatonic solids.

A neoplatonic solid is an undented polyhedron with equilateral triangle
faces, meeting at most six to a vertex.  It is prime if it has no
degree-3 vertex.

## Quick start

    make all VMAX=12       # everything through v=12 (~1 minute)

For large v, clone to a compute server:

    make all VMAX=50 JOBS=96

Requires: C compiler, python3, numpy, scipy.
GNU parallel is optional (used for faster sharded enumeration if present).

## Pipeline

    make seed-one  V=N      generate one seed file
    make seeds     VMAX=N   generate seed files through N
    make primes    VMAX=N   grow primes by recurrence: prime(v) = grow(prime(v-1)) ∪ seed(v)
    make solve     VMAX=N   Euclidean embeddings via homotopy from ideal
    make prove     VMAX=N   existence proofs (IEEE 754 + error bounds)
    make check     VMAX=N   independent validation (dent, length, defect, embed)
    make checkers           build the standalone check tools
    make all       VMAX=N   everything
    make status             show results

## Seed and prime semantics

- `seed(6)` = the octahedron (base case)
- `seed(v)` = empty for `v < 6` and `7 ≤ v < 12`
- `seed(v)` = fullerene duals via buckygen for `v ≥ 12`
- `prime(4)` = the tetrahedron (exceptional, does not seed the recurrence)
- `prime(5)` = empty
- `prime(v) = grow(prime(v-1)) ∪ seed(v)` for `v ≥ 6`

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
      prove.c                 existence prover (C, LAPACK SVD, IEEE 754 error bounds)
      prove_float.py          existence prover (Python, legacy)
      plantri_to_poly         planar_code → face-list converter
      dent_check.c            min link turning (undentedness)
      length_check.c          max edge length deviation
      defect_check.c          min angle defect
      embed_check.cpp         self-intersection test (requires CGAL, optional)
      check_all.c             runs dent + length + defect in one pass
    third_party/buckygen/   Brinkmann's fullerene generator (GPL)
    data/                   generated data (gitignored)
    run/                    scratch and logs (gitignored)
    PROOF_METHOD.txt        why the float prover is rigorous
