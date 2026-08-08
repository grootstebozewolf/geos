#!/usr/bin/env python3
# Adversarial hunt for GEOS CircularArcIntersector, modeled on
# NetTopologySuite.Proofs oracle/gen_arc_arc_tests.py + gen_arc_segment_tests.py.
# Exact-rational ground truth via fractions.Fraction; batched through one
# driver process for speed.
#
# Gated invariants (proven in the Rocq corpus; violations are real bugs):
#   ARC_ARC:     I1 on-both-circles, I2 shared-vertex>=1, I3 symmetry,
#                I4 reversal-stable, plus COINCIDENT only for identical circles
#   ARC_SEGMENT: I1 on-circle, I2 on-segment t in [0,1], I3 seg-reversal,
#                I4 arc-reversal, shared-vertex>=1 (noding robustness)
# Decisive exact-count mismatches (disc sign clear at float scale) are also
# reported; near-tangency (sub-ulp) divergence is only a signal.
import os
import random
import subprocess
import sys
from fractions import Fraction as F

DRIVER = os.environ.get("GEOS_DRIVER",
    os.path.join(os.path.dirname(__file__), "geos_arc_driver"))

random.seed(20260808)

# ---------- exact helpers (mirrors gen_arc_arc_tests.py) ----------

def circumcentre(p):
    (ax, ay), (bx, by), (cx, cy) = [(F(x), F(y)) for (x, y) in p]
    d = 2 * ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
    if d == 0:
        return None
    bk = bx * bx + by * by - ax * ax - ay * ay
    ck = cx * cx + cy * cy - ax * ax - ay * ay
    ox = ((cy - ay) * bk - (by - ay) * ck) / d
    oy = ((bx - ax) * ck - (cx - ax) * bk) / d
    r2 = (ox - ax) * (ox - ax) + (oy - ay) * (oy - ay)
    return (ox, oy, r2)

def exact_circle_count(a1, a2):
    o1, o2 = circumcentre(a1), circumcentre(a2)
    if o1 is None or o2 is None:
        return "DEGENERATE"
    (o1x, o1y, r1), (o2x, o2y, r2) = o1, o2
    dq = (o2x - o1x) ** 2 + (o2y - o1y) ** 2
    if dq == 0:
        return "COINCIDENT" if r1 == r2 else 0
    disc = 4 * dq * r1 - (dq + r1 - r2) ** 2
    return 2 if disc > 0 else (1 if disc == 0 else 0)

def exact_disc_rel(a1, a2):
    """Relative magnitude of the four-factor discriminant (0 if degenerate)."""
    o1, o2 = circumcentre(a1), circumcentre(a2)
    if o1 is None or o2 is None:
        return None
    (o1x, o1y, r1), (o2x, o2y, r2) = o1, o2
    dq = (o2x - o1x) ** 2 + (o2y - o1y) ** 2
    if dq == 0:
        return None
    disc = 4 * dq * r1 - (dq + r1 - r2) ** 2
    scale = (4 * dq * r1 + (dq + r1 - r2) ** 2) or F(1)
    return float(disc / scale)

def exact_line_count(arc, p, q):
    o = circumcentre(arc)
    if o is None:
        return "DEGENERATE"
    ox, oy, r2 = o
    (px, py), (qx, qy) = (F(p[0]), F(p[1])), (F(q[0]), F(q[1]))
    dx, dy = qx - px, qy - py
    len2 = dx * dx + dy * dy
    if len2 == 0:
        return "DEGENERATE"
    # squared distance from O to the supporting line
    num = (dy * (ox - px) - dx * (oy - py)) ** 2
    d2 = num / len2
    return 2 if r2 > d2 else (1 if r2 == d2 else 0)

def on_circle(pt, o, tol_rel=1e-9):
    ox, oy, r2 = float(o[0]), float(o[1]), float(o[2])
    d2 = (pt[0] - ox) ** 2 + (pt[1] - oy) ** 2
    return abs(d2 - r2) <= tol_rel * (1.0 + abs(r2)) * 64

# ---------- batched driver I/O ----------

class Batch:
    def __init__(self):
        self.stanzas = []
        self.meta = []

    def add_arc_arc(self, name, a1, a2):
        pts = list(a1) + list(a2)
        self.stanzas.append("ARC_ARC_XY\n" + "".join(f"{repr(x)} {repr(y)}\n" for x, y in pts))
        self.meta.append((name, ("AA", a1, a2)))

    def add_arc_seg(self, name, arc, p, q):
        pts = list(arc) + [p, q]
        self.stanzas.append("ARC_SEGMENT_XY\n" + "".join(f"{repr(x)} {repr(y)}\n" for x, y in pts))
        self.meta.append((name, ("AS", arc, p, q)))

    def run(self):
        proc = subprocess.run([DRIVER], input="".join(self.stanzas),
                              capture_output=True, text=True)
        lines = proc.stdout.strip().split("\n") if proc.stdout.strip() else []
        if len(lines) != len(self.meta):
            sys.exit(f"driver desync: {len(lines)} results for {len(self.meta)} cases\n"
                     + proc.stderr[-2000:])
        return list(zip(self.meta, lines))

def parse(line):
    tok = line.split()
    if not tok:
        return ("NAN", 0, [])
    if tok[0] in ("DEGENERATE", "COINCIDENT", "NAN"):
        return (tok[0], 0, [])
    n = int(tok[0])
    coords = [float.fromhex(t) if ("x" in t or "p" in t.lower()) else float(t)
              for t in tok[1:2 * n + 1]]
    return ("COUNT", n, [(coords[2 * i], coords[2 * i + 1]) for i in range(n)])

# ---------- case construction ----------

def rot_scale(pts, cs, sn, s, tx, ty):
    return tuple((s * (cs * x - sn * y) + tx, s * (sn * x + cs * y) + ty)
                 for x, y in pts)

def circle_point(cx, cy, r, theta):
    import math
    return (cx + r * math.cos(theta), cy + r * math.sin(theta))

violations = []
signals = []

def check_arc_arc(results):
    """First pass results: dict name->(kind,n,pts). Verify invariants."""
    by_name = {}
    for (name, spec), line in results:
        by_name[name] = (spec, parse(line))
    return by_name

# ============ build the hunt ============
batch = Batch()
cases = {}   # name -> spec for post-processing

def add_aa(name, a1, a2, shared=False):
    cases[name] = ("AA", a1, a2, shared)
    batch.add_arc_arc(name + "|base", a1, a2)
    batch.add_arc_arc(name + "|sym", a2, a1)
    batch.add_arc_arc(name + "|rev1", tuple(reversed(a1)), a2)
    batch.add_arc_arc(name + "|rev2", a1, tuple(reversed(a2)))

def add_as(name, arc, p, q, shared=False):
    cases[name] = ("AS", arc, p, q, shared)
    batch.add_arc_seg(name + "|base", arc, p, q)
    batch.add_arc_seg(name + "|segrev", arc, q, p)
    batch.add_arc_seg(name + "|arcrev", tuple(reversed(arc)), p, q)

import math

# --- A. a==0 family: d^2 + r1^2 == r2^2 (Pythagorean), radical line through c1
for (r1, d, r2) in [(3, 4, 5), (6, 8, 10), (5, 12, 13), (8, 15, 17), (20, 21, 29)]:
    a1 = ((r1, 0), (0, r1), (-r1, 0))
    a2 = ((d + r2, 0), (d, r2), (d - r2, 0))
    add_aa(f"pyth[{r1},{d},{r2}]-upper/upper", a1, a2)
    a1l = ((r1, 0), (0, -r1), (-r1, 0))
    add_aa(f"pyth[{r1},{d},{r2}]-lower/upper", a1l, a2)
    # swapped radii: a==0 relative to the *other* centre
    b1 = ((r2, 0), (0, r2), (-r2, 0))
    b2 = ((d + r1, 0), (d, r1), (d - r1, 0))
    add_aa(f"pyth-swap[{r2},{d},{r1}]-upper/upper", b1, b2)

# --- B. shared-vertex hunts (I2): circles through a common integer point,
# many geometries, all hinge combos.
shared_pts = [(5, 0), (0, 0), (-3, 4), (1e8, 1), (0.1, 0.7)]
for i in range(200):
    sp = random.choice(shared_pts)
    r1 = random.choice([1, 2, 3, 5, 7, 11, 0.5, 100, 1e6])
    r2 = random.choice([1, 2, 3, 5, 7, 11, 0.5, 100, 1e6])
    th1 = random.uniform(0, 2 * math.pi)
    th2 = random.uniform(0, 2 * math.pi)
    # centres placed so that each circle passes through sp exactly? Only
    # approximately: construct centre at distance r from sp via cos/sin, then
    # rebuild the arc THROUGH sp exactly by using sp as an endpoint.
    c1 = (sp[0] + r1 * math.cos(th1), sp[1] + r1 * math.sin(th1))
    c2 = (sp[0] + r2 * math.cos(th2), sp[1] + r2 * math.sin(th2))
    sw1 = random.uniform(0.3, 5.5)
    sw2 = random.uniform(0.3, 5.5)
    ang1 = math.atan2(sp[1] - c1[1], sp[0] - c1[0])
    ang2 = math.atan2(sp[1] - c2[1], sp[0] - c2[0])
    o1 = random.choice([1, -1])
    o2 = random.choice([1, -1])
    arc1 = (sp, circle_point(*c1, r1, ang1 + o1 * sw1 / 2), circle_point(*c1, r1, ang1 + o1 * sw1))
    arc2 = (sp, circle_point(*c2, r2, ang2 + o2 * sw2 / 2), circle_point(*c2, r2, ang2 + o2 * sw2))
    # hinge combos: start/start, start/end, end/start, end/end
    combo = random.randrange(4)
    if combo & 1:
        arc1 = tuple(reversed(arc1))
    if combo & 2:
        arc2 = tuple(reversed(arc2))
    add_aa(f"shared[{i}]", arc1, arc2, shared=True)

# --- C. near-tangency ladder (external + internal), decisive & sub-ulp
for k, off in [("m30", -2**-30), ("m45", -2**-45), ("0", 0.0),
               ("p45", 2**-45), ("p30", 2**-30), ("dec+", 1e-3), ("dec-", -1e-3)]:
    dx = 10.0 + off
    a1 = ((5, 0), (0, 5), (-5, 0))
    a2 = ((dx - 5, 0), (dx, 5), (dx + 5, 0))
    add_aa(f"ext-tan[{k}]", a1, a2)
    dxi = 2.0 + off
    b2 = ((dxi + 3, 0), (dxi, 3), (dxi - 3, 0))
    a1n = ((5, 0), (0, -5), (-5, 0))
    add_aa(f"int-tan[{k}]", a1, b2)
    add_aa(f"int-tan-low[{k}]", a1n, b2)

# --- D. random general-position arc pairs, exact-count comparison
for i in range(300):
    def rnd_arc():
        cx, cy = random.uniform(-20, 20), random.uniform(-20, 20)
        r = random.uniform(0.1, 15)
        t0 = random.uniform(0, 2 * math.pi)
        sw = random.uniform(0.2, 6.0) * random.choice([1, -1])
        return (circle_point(cx, cy, r, t0),
                circle_point(cx, cy, r, t0 + sw / 2),
                circle_point(cx, cy, r, t0 + sw))
    add_aa(f"rand[{i}]", rnd_arc(), rnd_arc())

# --- E. arc-segment: chords, tangents, shared endpoints, random
for i in range(200):
    cx, cy = random.uniform(-10, 10), random.uniform(-10, 10)
    r = random.choice([1, 3, 5, 0.25, 50])
    t0 = random.uniform(0, 2 * math.pi)
    sw = random.uniform(0.3, 6.0) * random.choice([1, -1])
    arc = (circle_point(cx, cy, r, t0),
           circle_point(cx, cy, r, t0 + sw / 2),
           circle_point(cx, cy, r, t0 + sw))
    kind = i % 4
    if kind == 0:      # segment with one endpoint exactly the arc start
        p = arc[0]
        q = (random.uniform(-15, 15), random.uniform(-15, 15))
        add_as(f"seg-shared[{i}]", arc, p, q, shared=True)
    elif kind == 1:    # secant through interior
        p = (cx - 2 * r, cy + random.uniform(-r / 2, r / 2))
        q = (cx + 2 * r, cy + random.uniform(-r / 2, r / 2))
        add_as(f"seg-secant[{i}]", arc, p, q)
    elif kind == 2:    # near-tangent horizontal line at y = cy + r +/- eps
        eps = random.choice([0.0, 2**-40, -2**-40, 1e-3, -1e-3])
        y = cy + r + eps
        add_as(f"seg-tan[{i}]", arc, (cx - 2 * r, y), (cx + 2 * r, y))
    else:              # fully random
        p = (random.uniform(-15, 15), random.uniform(-15, 15))
        q = (random.uniform(-15, 15), random.uniform(-15, 15))
        add_as(f"seg-rand[{i}]", arc, p, q)

# integer chord battery: exact tangent / secant / miss at integer scale
for r1 in (5, 25):
    arcU = ((r1, 0), (0, r1), (-r1, 0))
    add_as(f"chord-tan-int[{r1}]", arcU, (-2 * r1, r1), (2 * r1, r1))       # tangent at (0,r)
    add_as(f"chord-sec-int[{r1}]", arcU, (-2 * r1, 0), (2 * r1, 0))         # through (+-r,0)
    add_as(f"chord-miss-int[{r1}]", arcU, (-2 * r1, r1 + 1), (2 * r1, r1 + 1))

results = batch.run()
res = {name: parse(line) for (name, _spec), line in results}

# ============ evaluate ============
n_checked = 0
for name, spec in cases.items():
    n_checked += 1
    if spec[0] == "AA":
        _, a1, a2, shared = spec
        kind, n, pts = res[name + "|base"]
        ksym, nsym, _ = res[name + "|sym"]
        krev1, nrev1, _ = res[name + "|rev1"]
        krev2, nrev2, _ = res[name + "|rev2"]
        exact = exact_circle_count(a1, a2)
        rel = exact_disc_rel(a1, a2)

        # COINCIDENT must imply exactly-coincident circumcircles
        if kind == "COINCIDENT" and exact != "COINCIDENT":
            violations.append((name, f"BOGUS_COINCIDENT exact={exact} got={kind}"))
            continue
        if kind == "COUNT":
            o1, o2 = circumcentre(a1), circumcentre(a2)
            for pt in pts:
                if not (on_circle(pt, o1) and on_circle(pt, o2)):
                    violations.append((name, f"I1_OFF_CIRCLE pt={pt}"))
            if shared and n < 1:
                violations.append((name, f"I2_SHARED_VERTEX count={n}"))
            if ksym == "COUNT" and nsym != n:
                violations.append((name, f"I3_ASYMMETRY base={n} sym={nsym}"))
            if krev1 == "COUNT" and nrev1 != n:
                violations.append((name, f"I4_REV_ARC1 base={n} rev={nrev1}"))
            if krev2 == "COUNT" and nrev2 != n:
                violations.append((name, f"I4_REV_ARC2 base={n} rev={nrev2}"))
            # decisive count mismatch vs exact circle-pair count: only when
            # arcs cover the whole relevant region is count comparable; we
            # compare only 0-vs-nonzero when disc is decisively negative
            # (circles don't meet -> arcs can't meet).
            if isinstance(exact, int) and exact == 0 and n > 0 and rel is not None and rel < -1e-9:
                violations.append((name, f"PHANTOM_POINT exact_circles=0 got={n}"))
    else:
        _, arc, p, q, shared = spec
        kind, n, pts = res[name + "|base"]
        ksr, nsr, _ = res[name + "|segrev"]
        kar, nar, _ = res[name + "|arcrev"]
        exact = exact_line_count(arc, p, q)
        if kind == "COUNT":
            o = circumcentre(arc)
            for pt in pts:
                if o and not on_circle(pt, o):
                    violations.append((name, f"I1_OFF_CIRCLE pt={pt}"))
                # I2: on segment, t in [0,1] within tolerance
                dx, dy = q[0] - p[0], q[1] - p[1]
                L2 = dx * dx + dy * dy
                if L2 > 0:
                    t = ((pt[0] - p[0]) * dx + (pt[1] - p[1]) * dy) / L2
                    cross = (pt[0] - p[0]) * dy - (pt[1] - p[1]) * dx
                    if not (-1e-7 <= t <= 1 + 1e-7):
                        violations.append((name, f"I2_OFF_SEGMENT t={t}"))
                    if abs(cross) > 1e-6 * (1 + L2):
                        violations.append((name, f"I2_NOT_COLLINEAR cross={cross}"))
            if shared and n < 1:
                violations.append((name, f"SHARED_VERTEX_SEG count={n}"))
            if ksr == "COUNT" and nsr != n:
                violations.append((name, f"I3_SEGREV base={n} rev={nsr}"))
            if kar == "COUNT" and nar != n:
                violations.append((name, f"I4_ARCREV base={n} rev={nar}"))
            if isinstance(exact, int) and exact == 0 and n > 0:
                violations.append((name, f"PHANTOM_POINT exact_line=0 got={n}"))

print(f"# checked {n_checked} cases ({len(batch.meta)} driver invocations)")
if violations:
    print(f"# {len(violations)} VIOLATIONS:")
    for name, v in violations:
        base = res.get(name + "|base")
        print(f"!! [{name}] {v}   base_result={base}")
        spec = cases[name]
        if spec[0] == "AA":
            print(f"   a1={spec[1]} a2={spec[2]}")
        else:
            print(f"   arc={spec[1]} p={spec[2]} q={spec[3]}")
    sys.exit(1)
print("# all gated invariants hold")
