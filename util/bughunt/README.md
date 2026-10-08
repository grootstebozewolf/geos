# curve-bughunter: GEOS circular-arc code vs NetTopologySuite.Proofs

A differential tester for the GEOS curve code (`geom::CircularArc`,
`algorithm::CircularArcs`, `algorithm::CircularArcIntersector`,
`noding::GeometryNoder` on curves, and curve `getLength` / `getArea` /
`getLinearized`).  Random and edge-case inputs go through GEOS; every answer
is checked against an oracle that implements what NetTopologySuite.Proofs
proves, in exact rational arithmetic (Python `Fraction`) plus 60-digit
`mpmath` where a square root or trig function is unavoidable.

Nothing here changes GEOS. It is a test harness only.

## Build and run

```sh
# 1. build GEOS as usual (Release is fine)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j

# 2. build the probe (C++ driver that reads case lines on stdin, prints JSON)
util/bughunt/build.sh                      # GEOS_BUILD=... to point elsewhere

# 3. optional: build the NTS.Proofs extracted oracle (OCaml + zarith)
#    clones NTS.Proofs (pinned 57f82369) unless NTS_PROOFS points at a checkout;
#    needs extracted.ml from `make -f Makefile.gen` in that checkout
NTS_PROOFS=/path/to/NetTopologySuite.Proofs util/bughunt/build_oracle.sh

# 4. hunt (python3 with mpmath)
cd util/bughunt
python3 hunt.py run --seed 1 --max-cases 100000 --minutes 45 --out results
python3 hunt.py case 'AS 5 0 0 5 -5 0 3 0 3.000000001 10'   # rerun one case, verbose
python3 hunt.py shrink '<case line>'                       # minimise a failing case
```

`repros.txt` holds the minimal reproducer for every finding (see RESULTS.md);
`python3 replay.py results/findings.jsonl [check ...]` re-evaluates the cases of a
previous run, e.g. against a patched GEOS build, and reports how many still fail.

`--no-oracle` skips the OCaml oracle (the Python exact oracle is always used).
Output: `results/summary.json` (counts) and `results/findings.jsonl` (one line
per finding, with the case line that reproduces it).  The run is deterministic
for a given `--seed`.

Case line formats (see `geos_probe.cpp`):

| op | arguments | what is called |
|----|-----------|----------------|
| `A`  | `x0 y0 x1 y1 x2 y2 maxDev` | CircularArc orientation, centre, radius, angle, length, area, isCircle/isLinear, linearization (and the reversed arc) |
| `AA` | two arcs (12 numbers) | CircularArcIntersector both orders, distance, closestPoints |
| `AS` | arc + segment (10 numbers) | CircularArcIntersector arc/segment both segment directions, distance, closestPoints |
| `AP` | arc + point | distance, closestPoint, containsPoint |
| `SS` | two segments via the arc intersector | |
| `G`  | WKT | getLength, getArea, getLinearized |
| `N`  | `0|1` WKT | GeometryNoder (1 = preserve compound curves, which also sets onlyFirstGeomEdges as GeometrySplitter does) |

## What is checked, and which proof backs it

Checks marked **proved** compare against a property that has a Qed proof in
NetTopologySuite.Proofs (pin 57f82369; none of the files below contain
`Admitted`).  Checks marked **sanity only** are not backed by a proof; a
mismatch there is still reported, but labelled as such.

| check | property | backing |
|-------|----------|---------|
| centre, radius | circumcentre equidistant from the three controls | proved: `CircumcentreQSound.v` `circumcentre_formula_equidistant`; arc validity `CurveGeometry.v` |
| orientation | exact sign of the control triangle | proved: `ArcOrient.v` (b64-exact orientation) |
| sweep / length | length = r * sweep; full turn = 2*pi*r; chord <= length | proved: `ArcLength.v` `arc_length`, `chord_subtended_sq`, `chord_le_arc_length`, `arc_length_full_turn`; `ArcChordLength.v` `arc_chord_le_arc_length`; cos of the sweep: `AngleBetween.v` `cos_angle_between` (oracle side: `arc_invariants_q` in `oracle/driver.ml`) |
| segment area | (r^2/2)(theta - sin theta) >= 0, full disc = pi r^2 | proved: `ArcArea.v` `segment_area`, `segment_area_nonneg`, `segment_area_full_disc`; sign flip `CurvePolygonOrientation.v` `segment_area_neg` |
| curve-polygon area (`G` area) | assembly of segment areas into a polygon area | **sanity only** (the assembly is deferred in `CurvePolygonOrientation.v`) |
| linearization | vertex count n = max(2, ceil(sweep/step)); vertices on the circle, inside the span, monotone in angle; last vertex = p2 exactly; reverse commutes; chord sagitta <= maxDeviation | proved: `ArcLinearize.v` `subdiv_n_ok`, `lin_vertex_on_arc`, `lin_theta_mono_pos/neg`, `lin_full_circle_closed`, `lin_reverse_commutes`; `ArcLinearizeBound.v` `dev_to_step`, `chord_sagitta_le_dev`, `chord_approx_error_bound` |
| arc/arc points | radical-line points on both circles, kept iff inside both spans; shared endpoints always intersect; symmetric | proved: `ArcArcCircles.v` `radical_points_on_circles`, `two_circles_radical_point_unique`; `ArcSpanAtan2.v` `arc_span_contains_atan2_iff_chord_sign`, `arc_arc_intersects_of_atan2_radical_span`; `LECArcRow.v` `arc_span_contains_iff_sign`; `ArcArcSound.v` `arc_arc_intersects_shared_vertex`; `ArcIntersect.v` `arc_arc_intersects_sym` |
| arc/segment points | line-circle points, kept iff in the span and between P and Q | proved: `ArcSegmentCircles.v` `line_circle_radical_point`, `arc_line_circle_intersect` (the "between P Q" clamp is the definition in `ArcIntersect.v`) |
| point distance, containsPoint | arc_dist = min(foot if in span, endpoints) | proved: `LECArcRow.v` `arc_dist`, `arc_dist_exact`, `on_circle_side_zero_is_endpoint` |
| arc/arc, arc/segment distance | lower bounds d(O1,O2)-r1-r2 and d(O,line)-r | proved: `ArcArcDistance.v` `two_circles_dist_lower`, `arc_arc_dist_external`; `ArcSegmentDistance.v` `circle_line_dist_lower` |
| distance closest points lie on their arcs, distance <= endpoint pairs | definitional | **sanity only** |
| noding preserves length; every crossing is a node | split pieces cover the arc and are disjoint; length additive | proved: `ArcSplitAtNode.v` `arc_split_cover`, `arc_split_disjoint`; `CurveLength.v` `curve_length_additive` |
| full circle CIRCULARSTRING(A,B,A) | circle with diameter AB | convention, not a correctness theorem: ADR-0005 / `IntakeWalkerCircle.v`, parameterisation `IntakeCircle.v` `circle_full_param` |
| cocircular overlap (same circle) | overlap arcs lie on both inputs | **sanity only** |
| collinear "arcs" | not valid arcs in `CurveGeometry.v` | **sanity only** (internal consistency) |

## Tolerances and conditioning

Every float comparison uses `1e-9 * |value| + 2^-40 * M + 64 * S`, where M is
the largest coordinate magnitude and S is the measured sensitivity of the
exact answer to 2-ulp componentwise perturbations of the inputs.  Inputs
whose circle moves by more than `1e-9 * r` under such perturbations, near
tangencies whose gap is below `2^-36 * M`, and points within `1e-6 * r` of a
span boundary are classified `ill_conditioned` instead of as bugs.

The oracle_arc cross-check (`ORA.*`) compares the Python oracle with the
OCaml oracle built from NTS.Proofs' extraction; disagreements there are about
the oracle, not GEOS.
