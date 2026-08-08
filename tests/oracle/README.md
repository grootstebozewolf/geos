# Differential testing against the NetTopologySuite.Proofs oracle

This directory contains a differential-test harness that gates GEOS's
`CircularArcIntersector` against the formally verified reference oracle from
[NetTopologySuite.Proofs](https://github.com/grootstebozewolf/NetTopologySuite.Proofs)
(the RocqRefRunner `oracle_bin`, extracted from Rocq/Coq proofs) and against
exact rational arithmetic.

## Components

- `geos_arc_driver.cpp` — a GEOS-side driver speaking the oracle's
  stdin/stdout protocol for the `ARC_ARC_XY` and `ARC_SEGMENT_XY` modes.
  Because it speaks the same protocol, the Proofs repository's own gated
  generators (`oracle/gen_arc_arc_tests.py`, `oracle/gen_arc_segment_tests.py`)
  can run against GEOS directly via `ORACLE_BIN=<path-to-geos_arc_driver>`.
- `hunt.py` — adversarial invariant gate with exact-rational ground truth
  (`fractions.Fraction`), modeled on the Proofs generators. Gates the proven
  invariants: I1 on-circle(s), I2 shared-vertex ⇒ count ≥ 1, I3 symmetry,
  I4 reversal stability, plus exact-count comparison on decisive
  configurations.
- `hunt2.py` — focused shared-vertex robustness hunt: arcs tangent at a
  shared vertex (arc/arc) and segments leaving a shared endpoint
  tangentially (arc/segment), the classes where rounded early-exits used to
  discard exact shared vertices.
- `diff_oracle.py` — head-to-head comparison of GEOS vs `oracle_bin`
  (counts and coordinates, magnitude-scaled tolerance) on random and
  crossing-biased corpora.

## Building the driver

```sh
g++ -O2 -std=c++17 -I../../include -I../../build/include \
    geos_arc_driver.cpp -L../../build/lib -lgeos \
    -Wl,-rpath,$PWD/../../build/lib -o geos_arc_driver
```

## Running

```sh
python3 hunt.py                 # exact-rational invariant gate
python3 hunt2.py                # shared-vertex robustness gate
ORACLE_BIN=/path/to/oracle_bin python3 diff_oracle.py   # vs verified oracle
```

`oracle_bin` is published by the Proofs repository's `build-oracle.yml`
workflow as the `oracle-bin-linux` artifact, or can be built from source
(opam: `rocq-core.9.2.0` + `rocq-stdlib.9.2.0` + `coq-flocq.4.2.2`, then
`rocq makefile -f _CoqProject.full -o Makefile.gen &&
make -f Makefile.gen theories-flocq/Validate_binary64_extract.vo &&
make -C oracle`).

## Semantics notes

- The oracle reports `COINCIDENT` for arcs on exactly-coincident
  circumcircles without enumerating overlap arcs; GEOS enumerates the
  overlap (`COCIRCULAR_INTERSECTION`). The harness treats these as
  equivalent classifications.
- Near-tangency count divergence at sub-ulp discriminants is the documented
  float frontier of both implementations and is reported as a signal, not a
  failure.
- GEOS's exact shared-control-point pre-pass implements the proven
  `arc_arc_intersects_shared_vertex` invariant directly, and can therefore
  report an intersection at an exactly-shared vertex that the float oracle
  misses; `diff_oracle.py` corpora avoid bitwise-shared vertices for this
  reason (they are gated by `hunt2.py` instead).
