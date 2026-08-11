#!/usr/bin/env python3
# Head-to-head differential: GEOS CircularArcIntersector vs the verified
# RocqRefRunner oracle (NetTopologySuite.Proofs oracle_bin) on ARC_ARC_XY and
# ARC_SEGMENT_XY. Compares result kind, count, and coordinates (as sets,
# magnitude-scaled tolerance). Near-tangency count divergence (both counts
# plausible for sub-ulp discriminants) is reported as a signal, not a failure.
import math
import os
import random
import subprocess
import sys
from fractions import Fraction as F

GEOS = os.environ.get("GEOS_DRIVER",
    os.path.join(os.path.dirname(__file__), "geos_arc_driver"))
ORACLE = os.environ.get("ORACLE_BIN", "/workspace/nettopologysuite.proofs/oracle/oracle_bin")

random.seed(1729)

def circumcentre(p):
    (ax, ay), (bx, by), (cx, cy) = [(F(x), F(y)) for (x, y) in p]
    d = 2 * ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
    if d == 0:
        return None
    bk = bx * bx + by * by - ax * ax - ay * ay
    ck = cx * cx + cy * cy - ax * ax - ay * ay
    ox = ((cy - ay) * bk - (by - ay) * ck) / d
    oy = ((bx - ax) * ck - (cx - ax) * bk) / d
    r2 = (ox - ax) ** 2 + (oy - ay) ** 2
    return (ox, oy, r2)

def near_tangent_aa(a1, a2, rel=1e-12):
    o1, o2 = circumcentre(a1), circumcentre(a2)
    if o1 is None or o2 is None:
        return True
    (o1x, o1y, r1), (o2x, o2y, r2) = o1, o2
    dq = (o2x - o1x) ** 2 + (o2y - o1y) ** 2
    if dq == 0:
        return True
    disc = 4 * dq * r1 - (dq + r1 - r2) ** 2
    scale = 4 * dq * r1 + (dq + r1 - r2) ** 2
    return scale == 0 or abs(disc / scale) < rel

def near_tangent_as(arc, p, q, rel=1e-12):
    o = circumcentre(arc)
    if o is None:
        return True
    ox, oy, r2 = o
    (px, py), (qx, qy) = (F(p[0]), F(p[1])), (F(q[0]), F(q[1]))
    dx, dy = qx - px, qy - py
    len2 = dx * dx + dy * dy
    if len2 == 0:
        return True
    d2 = (dy * (ox - px) - dx * (oy - py)) ** 2 / len2
    scale = r2 + d2
    return scale == 0 or abs((r2 - d2) / scale) < rel

def circle_point(cx, cy, r, theta):
    return (cx + r * math.cos(theta), cy + r * math.sin(theta))

cases = []   # (name, kind, payload, stanza)

def add_aa(name, a1, a2):
    pts = list(a1) + list(a2)
    stanza = "ARC_ARC_XY\n" + "".join(f"{x!r} {y!r}\n" for x, y in pts)
    cases.append((name, "AA", (a1, a2), stanza))

def add_as(name, arc, p, q):
    pts = list(arc) + [p, q]
    stanza = "ARC_SEGMENT_XY\n" + "".join(f"{x!r} {y!r}\n" for x, y in pts)
    cases.append((name, "AS", (arc, p, q), stanza))

def rnd_arc(scale=20):
    cx, cy = random.uniform(-scale, scale), random.uniform(-scale, scale)
    r = random.uniform(0.05, scale)
    t0 = random.uniform(0, 2 * math.pi)
    sw = random.uniform(0.2, 6.0) * random.choice([1, -1])
    return (circle_point(cx, cy, r, t0),
            circle_point(cx, cy, r, t0 + sw / 2),
            circle_point(cx, cy, r, t0 + sw))

for i in range(3000):
    add_aa(f"aa[{i}]", rnd_arc(), rnd_arc())
for i in range(1500):
    # crossing-biased: two circles with overlapping extents
    cx, cy = random.uniform(-5, 5), random.uniform(-5, 5)
    r1 = random.uniform(1, 10)
    r2 = random.uniform(1, 10)
    d = random.uniform(abs(r1 - r2) * 0.9, (r1 + r2) * 1.1)
    th = random.uniform(0, 2 * math.pi)
    c2 = (cx + d * math.cos(th), cy + d * math.sin(th))
    def arc_on(c, r):
        t0 = random.uniform(0, 2 * math.pi)
        sw = random.uniform(0.5, 6.0) * random.choice([1, -1])
        return (circle_point(c[0], c[1], r, t0),
                circle_point(c[0], c[1], r, t0 + sw / 2),
                circle_point(c[0], c[1], r, t0 + sw))
    add_aa(f"aax[{i}]", arc_on((cx, cy), r1), arc_on(c2, r2))
for i in range(3000):
    arc = rnd_arc()
    p = (random.uniform(-25, 25), random.uniform(-25, 25))
    q = (random.uniform(-25, 25), random.uniform(-25, 25))
    add_as(f"as[{i}]", arc, p, q)

def run(bin_path):
    inp = "".join(st for (_, _, _, st) in cases)
    proc = subprocess.run([bin_path], input=inp, capture_output=True, text=True)
    lines = proc.stdout.strip().split("\n")
    if len(lines) != len(cases):
        sys.exit(f"{bin_path}: desync {len(lines)} vs {len(cases)}\n" + proc.stderr[-800:])
    return lines

def parse(line):
    tok = line.split()
    if not tok:
        return ("NAN", 0, [])
    if tok[0] in ("DEGENERATE", "NAN") or tok[0].startswith("COINCIDENT"):
        return (tok[0].split("=")[0], 0, [])
    n = int(tok[0])
    coords = [float.fromhex(t) if ("x" in t or "p" in t.lower()) else float(t)
              for t in tok[1:2 * n + 1]]
    return ("COUNT", n, [(coords[2 * i], coords[2 * i + 1]) for i in range(n)])

geos_lines = run(GEOS)
oracle_lines = run(ORACLE)

mismatch_count = 0
mismatch_coord = 0
signal = 0
worst = (0.0, None)
for (name, kind, payload, _), gl, ol in zip(cases, geos_lines, oracle_lines):
    gk, gn, gpts = parse(gl)
    ok_, on_, opts = parse(ol)
    if gk != ok_ or (gk == "COUNT" and gn != on_):
        nt = near_tangent_aa(*payload) if kind == "AA" else near_tangent_as(*payload)
        if nt:
            signal += 1
        else:
            mismatch_count += 1
            if mismatch_count <= 10:
                print(f"!! COUNT [{name}] geos={gl!r} oracle={ol!r}")
                print(f"   payload={payload}")
        continue
    if gk == "COUNT" and gn > 0:
        # compare as sets
        scale = max(1.0, *(abs(c) for pt in opts for c in pt))
        tol = 1e-9 * scale
        unmatched = list(opts)
        maxerr = 0.0
        for gp in gpts:
            best = min(unmatched, key=lambda op: (op[0]-gp[0])**2 + (op[1]-gp[1])**2)
            err = math.hypot(best[0]-gp[0], best[1]-gp[1])
            maxerr = max(maxerr, err)
            unmatched.remove(best)
        if maxerr > worst[0]:
            worst = (maxerr, name)
        if maxerr > tol:
            mismatch_coord += 1
            if mismatch_coord <= 10:
                print(f"!! COORD [{name}] err={maxerr:g} geos={gl!r} oracle={ol!r}")
                print(f"   payload={payload}")

print(f"# {len(cases)} cases: count_mismatch={mismatch_count} coord_mismatch={mismatch_coord} "
      f"near_tangency_signals={signal} worst_coord_err={worst[0]:g} ({worst[1]})")
sys.exit(1 if (mismatch_count or mismatch_coord) else 0)
