# Run of 2026-10-08 (Europe/Amsterdam), GEOS 0d790d1 (libgeos/geos main)

`python3 hunt.py run --seed 1 --max-cases 100000 --minutes 45` — 100,000 cases,
about 7 minutes. Counts are finding instances (one case can raise several).
Classes: `ill_conditioned` = the exact answer itself is unstable under 2-ulp input
perturbation (not counted as a bug); `ORA.*` = NTS.Proofs' own oracle_arc
disagrees with the exact oracle (an oracle issue, not GEOS).

| check | count | finding (see repros.txt) |
|-------|------:|---------|
| AS.missing / AS.extra / AS.inaccurate_point / AS.distance_nonzero / AS.closest_point_off | 1839 / 475 / 145 / 532 / 7 | 1 (circleIntersectsLine) and 2 (segment endpoint on arc) |
| N.missing_node | 219 | 1 |
| AS.segment_direction | 430 | 2 (and 1) |
| N.arc_not_split | 53 | 2 |
| AA.missing | 72 | 3 (T-junction), a few are 11 (precision) |
| AA.cocircular_missed | 335 | 4 |
| AA.cocircular_arc_off | 451 | 5 |
| AA.closest_point_off_arc | 202 | 6 |
| N.length_not_preserved | 565 | 7 (full circle) and 1 (wrong node positions) |
| N.sliver_piece | 134 | 2, 3 and 15 (unsnapped nodes) |
| AP.containsPoint_fullcircle | 201 | 8 |
| AS.closest_points_swapped | 227 | 9 |
| AA.count_type_mismatch / AS.count_type_mismatch | 94 / 102 | 10 |
| A.radius / A.sweep / A.length / A.segment_area / A.linear_for_valid_arc / A.chord_le_length | 1653 / 1750 / 1393 / 1128 (+707 precision) / 822 / 581 | 11 (near-collinear, large coordinates) and 12 (near-closed) |
| L.count / L.deviation / L.on_circle / L.span | 202 / 257 / 527 / 335 | 11 and 12 (linearization inherits the bad centre/sweep) |
| AP.distance, G.length | 119, 5 | 11 / 12 |
| N.exception | 105 | 14 (63), 16 (40), 15 (2) |
| ORA.ARC_ARC_XY / ORA.ARC_SEGMENT_XY | 1311 / 580 | oracle_arc misses exact shared vertices |
| ill_conditioned (all checks) | 3385 | not bugs |

Root causes were confirmed by temporarily patching GEOS (not committed) and
replaying the findings with `replay.py`:

* parametric, centre-relative line/circle intersection in `circleIntersectsLine`:
  AS.distance_nonzero 542 -> 0, N.missing_node 156 -> 0, AS.extra 403 -> 2,
  AS.missing 2227 -> 273;
* snapping computed points to segment endpoints within 1e-9 r: AS.segment_direction 429 -> 0,
  AS.missing 273 -> 184 (the remaining ones were ill-conditioned chords, now classified so);
* snapping arc/arc points to any arc endpoint: AA.missing 80 -> 11;
* treating |d| and |r1 - r2| below 1e-12 r as cocircular: AA.cocircular_missed 643 -> 50;
* `COUNTERCLOCKWISE` at CircularArcIntersector.cpp:407: AA.cocircular_arc_off fully fixed
  (with conditioning gating);
* wrapping theta >= 2 pi in `CircularArcs::containsAngle`: AA.closest_point_off_arc 221 -> 1;
* giving split pieces of a full circle CLOCKWISE orientation in `NodableArcString::getNoded`
  (together with the line/circle fix): N.length_not_preserved 552 -> 5.
