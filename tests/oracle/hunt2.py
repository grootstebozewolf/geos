#!/usr/bin/env python3
# Focused hunt: shared-vertex robustness against float early-exits.
#   F. arc-arc: circles externally/internally tangent AT the shared vertex;
#      float d vs r1+r2 (or |r1-r2|) can round outward -> early exit skips
#      the endpoint-equality snap -> count 0 despite exact shared vertex.
#   G. arc-segment: segment leaves the shared vertex radially (outward or
#      inward-grazing) -> circle-line discriminant ~ 0 +/- noise; a negative
#      rounding skips the endpoint snap entirely.
import math
import os
import random
import subprocess
import sys

DRIVER = os.environ.get("GEOS_DRIVER",
    os.path.join(os.path.dirname(__file__), "geos_arc_driver"))

random.seed(424242)

def circle_point(cx, cy, r, theta):
    return (cx + r * math.cos(theta), cy + r * math.sin(theta))

stanzas = []
meta = []

def add_aa(name, a1, a2):
    pts = list(a1) + list(a2)
    stanzas.append("ARC_ARC_XY\n" + "".join(f"{x!r} {y!r}\n" for x, y in pts))
    meta.append((name, "AA", (a1, a2)))

def add_as(name, arc, p, q):
    pts = list(arc) + [p, q]
    stanzas.append("ARC_SEGMENT_XY\n" + "".join(f"{x!r} {y!r}\n" for x, y in pts))
    meta.append((name, "AS", (arc, p, q)))

# --- F. tangent-at-shared-vertex arc pairs, random geometry ---
for i in range(4000):
    # tangency point V; circles tangent at V along random direction u.
    vx, vy = random.choice([(0, 0), (5, 0), (0.1, 0.7), (1e7, 3), (-3, 4)])
    th = random.uniform(0, 2 * math.pi)
    r1 = random.uniform(0.01, 1000)
    ext = random.random() < 0.5
    r2 = random.uniform(0.01, 1000)
    ux, uy = math.cos(th), math.sin(th)
    c1 = (vx - r1 * ux, vy - r1 * uy)
    # external: centres on opposite sides; internal: same side
    c2 = (vx + r2 * ux, vy + r2 * uy) if ext else (vx - r2 * ux, vy - r2 * uy)
    if not ext and abs(r1 - r2) < 1e-3:
        continue
    a1ang = math.atan2(vy - c1[1], vx - c1[0])
    a2ang = math.atan2(vy - c2[1], vx - c2[0])
    sw1 = random.uniform(0.3, 4.0) * random.choice([1, -1])
    sw2 = random.uniform(0.3, 4.0) * random.choice([1, -1])
    arc1 = ((vx, vy),
            circle_point(*c1, r1, a1ang + sw1 / 2),
            circle_point(*c1, r1, a1ang + sw1))
    arc2 = ((vx, vy),
            circle_point(*c2, r2, a2ang + sw2 / 2),
            circle_point(*c2, r2, a2ang + sw2))
    if i % 2:
        arc1 = tuple(reversed(arc1))
    add_aa(f"F[{i}]{'ext' if ext else 'int'}", arc1, arc2)

# --- G. segment leaving shared vertex radially ---
for i in range(4000):
    cx, cy = random.uniform(-10, 10), random.uniform(-10, 10)
    r = random.uniform(0.01, 1000)
    t0 = random.uniform(0, 2 * math.pi)
    sw = random.uniform(0.3, 4.0) * random.choice([1, -1])
    arc = (circle_point(cx, cy, r, t0),
           circle_point(cx, cy, r, t0 + sw / 2),
           circle_point(cx, cy, r, t0 + sw))
    v = arc[0] if i % 2 == 0 else arc[2]
    # direction: exactly radial outward, or tangent, or slight graze
    mode = i % 3
    rx, ry = v[0] - cx, v[1] - cy
    L = math.hypot(rx, ry)
    if mode == 0:      # radially outward
        d = (rx / L, ry / L)
    elif mode == 1:    # tangent at v
        d = (-ry / L, rx / L)
    else:              # tangent + tiny outward tilt
        eps = 2.0 ** random.randint(-45, -20)
        d = (-ry / L + eps * rx / L, rx / L + eps * ry / L)
    ln = random.uniform(0.1, 10)
    q = (v[0] + ln * d[0], v[1] + ln * d[1])
    add_as(f"G[{i}]m{mode}", arc, v, q)

proc = subprocess.run([DRIVER], input="".join(stanzas), capture_output=True, text=True)
lines = proc.stdout.strip().split("\n")
if len(lines) != len(meta):
    sys.exit(f"desync {len(lines)} vs {len(meta)}: {proc.stderr[-500:]}")

fails = []
for (name, kind, spec), line in zip(meta, lines):
    tok = line.split()
    if tok[0] in ("DEGENERATE", "NAN"):
        continue
    if tok[0].startswith("COINCIDENT"):
        continue
    n = int(tok[0])
    if n < 1:
        fails.append((name, kind, spec, line))

print(f"# {len(meta)} cases, {len(fails)} shared-vertex misses")
for name, kind, spec, line in fails[:15]:
    print(f"!! [{name}] -> {line}")
    print(f"   spec={spec}")
if fails:
    sys.exit(1)
