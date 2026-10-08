#!/usr/bin/env python3
"""Differential bug hunter: GEOS circular-curve code vs NetTopologySuite.Proofs.

  python3 hunt.py run   --seed 1 --minutes 40 [--max-cases N] [--out results]
  python3 hunt.py case  'AA 0 0 5 5 10 0 0 5 5 0 10 5'      # rerun one case, verbose
  python3 hunt.py shrink 'A 0 1 0 -1 1e-17 1 0.1'           # minimise a failing case

Case lines are exactly what bin/geos_probe reads (see geos_probe.cpp).
"""
import argparse, collections, json, math, os, random, subprocess, sys, time
from fractions import Fraction as Fr
import mpmath as mp
import nts_oracle as O

HERE = os.path.dirname(os.path.abspath(__file__))
PROBE = os.path.join(HERE, 'bin', 'geos_probe')
ORACLE = os.path.join(HERE, 'bin', 'oracle_arc')
TWO_PI = 2 * math.pi

# ---------------------------------------------------------------------------
# processes
# ---------------------------------------------------------------------------

def run_probe(lines):
    p = subprocess.run([PROBE], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    out = p.stdout.splitlines()
    if p.returncode != 0 or len(out) != len(lines):
        # a crash: bisect to find the case(s) that kill the probe
        res = []
        for ln in lines:
            r = subprocess.run([PROBE], input=ln + '\n', capture_output=True, text=True)
            if r.returncode != 0 or not r.stdout.strip():
                res.append({'op': 'CRASH', 'rc': r.returncode, 'stderr': r.stderr[-300:]})
            else:
                res.append(json.loads(r.stdout.splitlines()[0]))
        return res
    return [json.loads(x) for x in out]


def fx(v):
    return float.hex(float(v))


def run_oracle(queries):
    """queries: list of (mode, [points]) -> list of output lines."""
    buf = []
    for mode, pts in queries:
        buf.append(mode)
        buf.extend('%s %s' % (fx(p[0]), fx(p[1])) for p in pts)
    p = subprocess.run([ORACLE], input='\n'.join(buf) + '\n', capture_output=True, text=True)
    out = p.stdout.splitlines()
    if len(out) != len(queries):
        raise RuntimeError('oracle_arc output mismatch: %d vs %d' % (len(out), len(queries)))
    return out


def parse_hex(tok):
    return float.fromhex(tok)


def fnum(v):
    if isinstance(v, str):
        return float(v)
    return float(v)

# ---------------------------------------------------------------------------
# case encoding
# ---------------------------------------------------------------------------

def r17(x):
    return repr(float(x))


def line_A(a, dev):
    return 'A ' + ' '.join(r17(c) for p in a for c in p) + ' ' + r17(dev)


def line_AA(a, b):
    return 'AA ' + ' '.join(r17(c) for p in a + b for c in p)


def line_AS(a, s):
    return 'AS ' + ' '.join(r17(c) for p in a + s for c in p)


def line_AP(a, p):
    return 'AP ' + ' '.join(r17(c) for c in [*a[0], *a[1], *a[2], *p])


def parse_case(line):
    t = line.split(None, 1)
    op = t[0]
    if op in ('G', 'N'):
        return op, line
    v = [float(x) for x in t[1].split()]
    P = [(v[i], v[i + 1]) for i in range(0, len(v) - (len(v) % 2), 2)]
    if op == 'A':
        return op, (P[0:3], v[6])
    if op == 'AA':
        return op, (P[0:3], P[3:6])
    if op == 'AS':
        return op, (P[0:3], P[3:5])
    if op == 'AP':
        return op, (P[0:3], P[3])
    if op == 'SS':
        return op, (P[0:2], P[2:4])
    raise ValueError(op)

# ---------------------------------------------------------------------------
# generators
# ---------------------------------------------------------------------------

PYTH_R = [5, 10, 13, 25, 50, 65, 85, 125, 325, 1105]
_pyth = {}


def pyth_points(r):
    if r not in _pyth:
        pts = set()
        for x in range(0, r + 1):
            y2 = r * r - x * x
            y = math.isqrt(y2)
            if y * y == y2:
                for sx in (1, -1):
                    for sy in (1, -1):
                        pts.add((sx * x, sy * y))
        _pyth[r] = sorted(pts)
    return _pyth[r]


def circle_pt(c, r, th):
    return (c[0] + r * math.cos(th), c[1] + r * math.sin(th))


def ulp_nudge(x, k):
    for _ in range(abs(k)):
        x = math.nextafter(x, math.inf if k > 0 else -math.inf)
    return x


class Gen:
    def __init__(self, rng):
        self.R = rng

    def coord_scale(self):
        return self.R.choice([1, 10, 1000, 1e6])

    def offset(self):
        return self.R.choice([0, 0, 0, 1e3, 1e6, 1e9]) * self.R.choice([1, -1])

    def rand_arc(self):
        R = self.R
        mode = R.random()
        if mode < 0.25:  # integer coordinates
            S = R.choice([10, 1000, 10 ** 6])
            while True:
                a = [(R.randint(-S, S), R.randint(-S, S)) for _ in range(3)]
                if O.cross(*[O.q(p) for p in a]) != 0:
                    return [tuple(map(float, p)) for p in a], 'int'
        if mode < 0.55:  # float arc from centre / radius / angles
            c = (R.uniform(-1, 1) * self.coord_scale() + self.offset(), R.uniform(-1, 1) * self.coord_scale())
            r = 10 ** R.uniform(-6, 8)
            t0 = R.uniform(-math.pi, math.pi)
            sw = R.choice([1, -1]) * R.uniform(1e-6, TWO_PI - 1e-6)
            tm = t0 + sw * R.uniform(0.05, 0.95)
            return [circle_pt(c, r, t0), circle_pt(c, r, tm), circle_pt(c, r, t0 + sw)], 'float'
        if mode < 0.65:  # Pythagorean (exactly concyclic integer points)
            r = R.choice(PYTH_R)
            pts = R.sample(pyth_points(r), 3)
            c = (R.randint(-50, 50), R.randint(-50, 50))
            return [(float(p[0] + c[0]), float(p[1] + c[1])) for p in pts], 'pyth'
        if mode < 0.72:  # half / quarter circles
            r = R.choice([1, 3, 5, 10, 1e3]); c = (R.choice([0, 7, -3, 1e6]), R.choice([0, 2, 1e6]))
            k = R.choice(['half', 'quarter', 'threequarter'])
            s = R.choice([1, -1])
            if k == 'half':
                return [(c[0] - r, c[1]), (c[0], c[1] + s * r), (c[0] + r, c[1])], 'half'
            if k == 'quarter':
                return [(c[0] + r, c[1]), circle_pt(c, r, s * math.pi / 4), (c[0], c[1] + s * r)], 'quarter'
            return [(c[0] + r, c[1]), (c[0] - r, c[1]), (c[0], c[1] - s * r)], 'threequarter'
        if mode < 0.80:  # near-collinear
            p0 = (R.uniform(-1, 1) * 100, R.uniform(-1, 1) * 100)
            p2 = (R.uniform(-1, 1) * 100, R.uniform(-1, 1) * 100)
            t = R.choice([0.5, R.uniform(0.01, 0.99), R.uniform(1.1, 3), R.uniform(-2, -0.1)])
            eps = 10 ** R.uniform(-16, -6) * R.choice([1, -1])
            dx, dy = p2[0] - p0[0], p2[1] - p0[1]
            p1 = (p0[0] + t * dx - eps * dy, p0[1] + t * dy + eps * dx)
            return [p0, p1, p2], 'nearcollinear'
        if mode < 0.90:  # p0 ~ p2 (near full circle) and exact full circles
            c = (R.uniform(-1, 1) * self.coord_scale(), R.uniform(-1, 1) * self.coord_scale())
            r = 10 ** R.uniform(-3, 6)
            t0 = R.uniform(-math.pi, math.pi)
            p0 = circle_pt(c, r, t0)
            p1 = circle_pt(c, r, t0 + R.choice([1, -1]) * R.uniform(0.5, 5.5))
            if R.random() < 0.4:
                return [p0, p1, p0], 'fullcircle'
            k = R.choice([1, 2, 3, 8, 64, 1024])
            p2 = (ulp_nudge(p0[0], R.choice([-k, 0, k])), ulp_nudge(p0[1], R.choice([-k, k])))
            return [p0, p1, p2], 'nearfull'
        # p1 near p0
        c = (R.uniform(-10, 10), R.uniform(-10, 10)); r = 10 ** R.uniform(-2, 4)
        t0 = R.uniform(-math.pi, math.pi); sw = R.choice([1, -1]) * R.uniform(0.1, 6)
        tm = t0 + sw * 10 ** R.uniform(-12, -4)
        return [circle_pt(c, r, t0), circle_pt(c, r, tm), circle_pt(c, r, t0 + sw)], 'p1nearp0'

    def arc_on_circle(self, c, r, exact=False):
        R = self.R
        if exact:
            pts = R.sample(pyth_points(r), 3)
            return [(float(p[0] + c[0]), float(p[1] + c[1])) for p in pts]
        t0 = R.uniform(-math.pi, math.pi); sw = R.choice([1, -1]) * R.uniform(0.05, TWO_PI - 0.05)
        return [circle_pt(c, r, t0), circle_pt(c, r, t0 + sw * R.uniform(0.1, 0.9)), circle_pt(c, r, t0 + sw)]

    def pair(self):
        R = self.R
        m = R.random()
        if m < 0.3:
            a, _ = self.rand_arc(); b, _ = self.rand_arc()
            # pull b near a so that they interact
            ca = O.circumcentre(*[O.q(p) for p in a])
            if ca is not None and R.random() < 0.8:
                cb = O.circumcentre(*[O.q(p) for p in b])
                if cb is not None:
                    sh = (float(ca[0] - cb[0]) + R.uniform(-1, 1) * abs(a[0][0] - float(ca[0])),
                          float(ca[1] - cb[1]) + R.uniform(-1, 1) * abs(a[0][1] - float(ca[1])))
                    b = [(p[0] + sh[0], p[1] + sh[1]) for p in b]
            return a, b, 'random'
        if m < 0.45:  # exact Pythagorean circles: crossings at integer points, tangency, coincident
            r1 = R.choice(PYTH_R[:7]); r2 = R.choice(PYTH_R[:7])
            c1 = (0, 0)
            k = R.choice(['tan_ext', 'tan_int', 'cross', 'coincident', 'concentric', 'shared'])
            if k == 'tan_ext':
                c2 = (r1 + r2, 0)
            elif k == 'tan_int':
                c2 = (abs(r1 - r2), 0) if r1 != r2 else (0, 0)
            elif k == 'coincident':
                r2 = r1; c2 = c1
            elif k == 'concentric':
                c2 = c1
            elif k == 'shared':
                p = R.choice(pyth_points(r1)); q2 = R.choice(pyth_points(r2))
                c2 = (p[0] - q2[0], p[1] - q2[1])
            else:
                c2 = (R.randint(-r1 - r2, r1 + r2), R.randint(-r1 - r2, r1 + r2))
            sh = (R.choice([0, 0, 17, 1e6]), R.choice([0, 0, -3, 1e6]))
            a = self.arc_on_circle(c1, r1, exact=True)
            b = self.arc_on_circle(c2, r2, exact=True)
            if k == 'shared':
                b = [(float(p[0] + c2[0]), float(p[1] + c2[1])) for p in R.sample(pyth_points(r2), 2)]
                b = [b[0], b[1], (float(p[0]), float(p[1]))] if R.random() < 0.5 else [(float(p[0]), float(p[1])), b[0], b[1]]
            if R.random() < 0.3:  # full circle variant
                a = [a[0], (2 * c1[0] - a[0][0], 2 * c1[1] - a[0][1]), a[0]]
            a = [(x + sh[0], y + sh[1]) for x, y in a]
            b = [(x + sh[0], y + sh[1]) for x, y in b]
            return a, b, 'pyth_' + k
        if m < 0.7:  # float tangent / near tangent / concentric / coincident
            c1 = (R.uniform(-5, 5) * self.coord_scale(), R.uniform(-5, 5)); r1 = 10 ** R.uniform(-3, 5)
            r2 = r1 * 10 ** R.uniform(-2, 2)
            th = R.uniform(-math.pi, math.pi)
            k = R.choice(['tan_ext', 'tan_int', 'near_tan', 'coincident', 'concentric', 'cross_end'])
            if k == 'tan_ext':
                d = r1 + r2
            elif k == 'tan_int':
                d = abs(r1 - r2)
            elif k == 'near_tan':
                d = (r1 + r2) * (1 + R.choice([1, -1]) * 10 ** R.uniform(-15, -6))
            elif k == 'coincident':
                d, r2 = 0.0, r1
            elif k == 'concentric':
                d = 0.0
            else:
                d = R.uniform(abs(r1 - r2), r1 + r2)
            c2 = (c1[0] + d * math.cos(th), c1[1] + d * math.sin(th))
            a = self.arc_on_circle(c1, r1)
            b = self.arc_on_circle(c2, r2)
            if k in ('tan_ext', 'tan_int', 'near_tan'):
                # make both arcs cover the tangent direction
                t1 = th if (k != 'tan_int' or r1 >= r2) else th + math.pi
                a = [circle_pt(c1, r1, t1 - 0.5), circle_pt(c1, r1, t1 + 0.1), circle_pt(c1, r1, t1 + 0.6)]
                t2 = th + math.pi if k != 'tan_int' else (th if r1 >= r2 else th + math.pi)
                if k == 'tan_int' and r1 >= r2:
                    t2 = th
                b = [circle_pt(c2, r2, t2 - 0.4), circle_pt(c2, r2, t2 - 0.1), circle_pt(c2, r2, t2 + 0.7)]
            if k == 'cross_end':
                # b starts exactly at a's end point, which is generally not on b's circle exactly
                b = [a[2], b[1], b[2]]
            if k == 'coincident' and R.random() < 0.5:
                # share an endpoint exactly
                b = [a[2], b[1], b[2]]
            return a, b, 'f_' + k
        if m < 0.85:  # arcs crossing at a shared endpoint / control point
            a, _ = self.rand_arc()
            b, _ = self.rand_arc()
            j = R.choice([0, 2])
            b = [a[j], b[1], b[2]] if R.random() < 0.5 else [b[0], b[1], a[j]]
            return a, b, 'shared_vertex'
        # cocircular overlapping arcs on the same float circle (same 3-pt circle hard -> reuse points)
        a, _ = self.rand_arc()
        if R.random() < 0.5:
            b = [a[1], a[2], a[0]]  # same circle exactly (same three points)
        else:
            b = [a[2], a[1], a[0]]  # reversed
        return a, b, 'cocircular'

    def seg_for(self, a, kind=None):
        R = self.R
        c = O.circumcentre(*[O.q(p) for p in a])
        if c is None:
            c = (Fr(a[1][0]), Fr(a[1][1]))
        c = (float(c[0]), float(c[1]))
        r = math.hypot(a[0][0] - c[0], a[0][1] - c[1]) or 1.0
        k = kind or R.choice(['through', 'endpoint_on_arc', 'tangent', 'vertical', 'near_vertical', 'chord', 'random', 'ctrl'])
        if k == 'through':
            p = circle_pt(c, r * R.uniform(0, 0.99), R.uniform(-math.pi, math.pi))
            th = R.uniform(-math.pi, math.pi); L = r * R.uniform(0.5, 4)
            return [(p[0] - L * math.cos(th), p[1] - L * math.sin(th)), (p[0] + L * math.cos(th), p[1] + L * math.sin(th))], k
        if k == 'endpoint_on_arc':
            e = R.choice([a[0], a[2]])
            th = R.uniform(-math.pi, math.pi); L = r * R.uniform(0.1, 3)
            other = (e[0] + L * math.cos(th), e[1] + L * math.sin(th))
            return ([e, other] if R.random() < 0.5 else [other, e]), k
        if k == 'ctrl':
            e = a[1]
            th = R.uniform(-math.pi, math.pi); L = r * R.uniform(0.1, 3)
            return [(e[0] - L * math.cos(th), e[1] - L * math.sin(th)), (e[0] + L * math.cos(th), e[1] + L * math.sin(th))], k
        if k == 'tangent':
            th = R.uniform(-math.pi, math.pi)
            t = circle_pt(c, r * (1 + R.choice([0, 0, 1e-12, -1e-12, 1e-8])), th)
            L = r * R.uniform(0.2, 2)
            return [(t[0] - L * math.sin(th), t[1] + L * math.cos(th)), (t[0] + L * math.sin(th), t[1] - L * math.cos(th))], k
        if k in ('vertical', 'near_vertical'):
            x = c[0] + r * R.uniform(-0.999, 0.999)
            dxv = 0.0 if k == 'vertical' else r * 10 ** R.uniform(-14, -6)
            return [(x, c[1] - 2 * r), (x + dxv, c[1] + 2 * r)], k
        if k == 'chord':
            return [a[0], a[2]], k
        P = (c[0] + R.uniform(-2, 2) * r, c[1] + R.uniform(-2, 2) * r)
        Q = (c[0] + R.uniform(-2, 2) * r, c[1] + R.uniform(-2, 2) * r)
        return [P, Q], k

    def pyth_seg(self):
        """Integer arc + integer segment ending exactly on the arc."""
        R = self.R
        r = R.choice(PYTH_R[:8]); pts = pyth_points(r)
        a = [tuple(map(float, p)) for p in R.sample(pts, 3)]
        e = R.choice(pts)
        o = (R.randint(-2 * r, 2 * r), R.randint(-2 * r, 2 * r))
        s = [tuple(map(float, e)), tuple(map(float, o))]
        if R.random() < 0.5:
            s.reverse()
        if R.random() < 0.3:  # segment through two concyclic points
            e2 = R.choice(pts); s = [tuple(map(float, e)), tuple(map(float, e2))]
        return a, s

    def curve_wkt(self):
        """Random CIRCULARSTRING / COMPOUNDCURVE / LINESTRING plus the arcs it holds."""
        R = self.R
        n = R.choice([1, 1, 2, 3])
        a, _ = self.rand_arc()
        if R.random() < 0.5:
            r = R.choice(PYTH_R[:6]); a = [tuple(map(float, p)) for p in R.sample(pyth_points(r), 3)]
        pts = list(a)
        for _ in range(n - 1):
            nxt, _ = self.rand_arc()
            sh = (pts[-1][0] - nxt[0][0], pts[-1][1] - nxt[0][1])
            pts += [(p[0] + sh[0], p[1] + sh[1]) for p in nxt[1:]]
        return pts

    def geometry(self):
        R = self.R
        k = R.choice(['cs', 'cc', 'cp', 'fullcircle_cp', 'cs_closed'])
        if k == 'cs':
            pts = self.curve_wkt()
            return 'G CIRCULARSTRING(%s)' % ', '.join('%r %r' % p for p in pts), k
        if k == 'cc':
            pts = self.curve_wkt()
            tail = (pts[-1][0] + R.uniform(-5, 5), pts[-1][1] + R.uniform(-5, 5))
            return 'G COMPOUNDCURVE(CIRCULARSTRING(%s), (%r %r, %r %r))' % (
                ', '.join('%r %r' % p for p in pts), pts[-1][0], pts[-1][1], tail[0], tail[1]), k
        if k == 'fullcircle_cp':
            r = R.choice([1, 2.5, 10, 1e3, 1e-3]); c = (R.uniform(-100, 100), R.uniform(-100, 100))
            t = R.uniform(-math.pi, math.pi)
            p0 = circle_pt(c, r, t); p1 = (2 * c[0] - p0[0], 2 * c[1] - p0[1])
            return 'G CURVEPOLYGON(CIRCULARSTRING(%r %r, %r %r, %r %r))' % (*p0, *p1, *p0), k
        # closed ring made of two arcs (lens / circle split in two)
        r = R.choice(PYTH_R[:6]); pts = R.sample(pyth_points(r), 4)
        pts.sort(key=lambda p: math.atan2(p[1], p[0]))
        ring = [pts[0], pts[1], pts[2], pts[3], pts[0]]
        ring = [tuple(map(float, p)) for p in ring]
        if k == 'cs_closed':
            return 'G CIRCULARSTRING(%s)' % ', '.join('%r %r' % p for p in ring), k
        return 'G CURVEPOLYGON(CIRCULARSTRING(%s))' % ', '.join('%r %r' % p for p in ring), k

    def noding(self):
        R = self.R
        if R.random() < 0.08:
            # full circle CIRCULARSTRING(A,B,A) touched by a line at 1 or 2 points (Pythagorean => exact)
            r = R.choice(PYTH_R[:6]); P = [tuple(map(float, p)) for p in pyth_points(r)]
            A = R.choice(P); B = (-A[0], -A[1]); T = R.choice(P)
            if R.random() < 0.5:
                ln = 'LINESTRING(%r %r, %r %r)' % (T[0], T[1], T[0] * 2, T[1] * 2)
            else:
                U = R.choice(P)
                ln = 'LINESTRING(%r %r, %r %r)' % (T[0], T[1], U[0], U[1])
            return 'N %d GEOMETRYCOLLECTION(CIRCULARSTRING(%r %r, %r %r, %r %r), %s)' % (
                R.choice([0, 1]), A[0], A[1], B[0], B[1], A[0], A[1], ln), 'node_fullcircle_touch'
        a, b, tag = self.pair()
        def wkt(arc):
            if R.random() < 0.2:
                return 'COMPOUNDCURVE(CIRCULARSTRING(%s), (%r %r, %r %r))' % (
                    ', '.join('%r %r' % p for p in arc), arc[2][0], arc[2][1], arc[2][0] + 1.5, arc[2][1] - 2.25)
            return 'CIRCULARSTRING(%s)' % ', '.join('%r %r' % p for p in arc)
        g2 = wkt(b)
        if R.random() < 0.3:
            s, _ = self.seg_for(a)
            g2 = 'LINESTRING(%r %r, %r %r)' % (*s[0], *s[1])
        pres = R.choice([0, 1])
        return 'N %d GEOMETRYCOLLECTION(%s, %s)' % (pres, wkt(a), g2), 'node_' + tag

    def case(self):
        R = self.R
        m = R.random()
        if m < 0.30:
            a, tag = self.rand_arc()
            ex = O.ExactArc(*a)
            rr = float(ex.r) if ex.kind != 'degen' else math.hypot(a[0][0] - a[2][0], a[0][1] - a[2][1])
            dev = R.choice([1e-4, 1e-3, 1e-2, 0.1, 0.5]) * (rr or 1.0)
            return line_A(a, dev), 'A_' + tag
        if m < 0.55:
            a, b, tag = self.pair()
            return line_AA(a, b), 'AA_' + tag
        if m < 0.75:
            if R.random() < 0.3:
                a, s = self.pyth_seg(); tag = 'pyth_end'
            else:
                a, _ = self.rand_arc(); s, tag = self.seg_for(a)
            return line_AS(a, s), 'AS_' + tag
        if m < 0.85:
            a, _ = self.rand_arc()
            k = R.random()
            if k < 0.3:
                r = R.choice(PYTH_R[:8]); pts = pyth_points(r)
                a = [tuple(map(float, p)) for p in R.sample(pts, 3)]
                p = tuple(map(float, R.choice(pts))); tag = 'on_circle'
            elif k < 0.4:
                c = O.circumcentre(*[O.q(x) for x in a]); p = (float(c[0]), float(c[1])) if c else a[1]; tag = 'centre'
            else:
                s = max(abs(c) for x in a for c in x) or 1
                p = (a[1][0] + R.uniform(-2, 2) * s, a[1][1] + R.uniform(-2, 2) * s); tag = 'random'
            return line_AP(a, p), 'AP_' + tag
        if m < 0.92:
            return self.geometry()
        return self.noding()

# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------

Finding = collections.namedtuple('Finding', 'cls check msg')


def mpd(p, q):
    return mp.sqrt((mp.mpf(p[0]) - q[0]) ** 2 + (mp.mpf(p[1]) - q[1]) ** 2)


def perturbed_sensitivity(pts, fn, rng, k=3):
    """Exact-quantity change under ~1-ulp input perturbations (conditioning)."""
    base = fn(pts)
    if base is None:
        return mp.inf
    worst = mp.mpf(0)
    for _ in range(k):
        # componentwise relative perturbation (a few ulps of each coordinate)
        q = [(x * (1 + rng.choice([-1, 1]) * 2.0 ** -51), y * (1 + rng.choice([-1, 1]) * 2.0 ** -51)) for x, y in pts]
        try:
            v = fn(q)
        except Exception:
            return mp.inf
        if v is None:
            return mp.inf
        worst = max(worst, abs(v - base))
    return worst


_kcache = {}


def arc_kappa(pts, rng):
    """Relative circle conditioning: how far (centre, radius) move, relative to r,
    under componentwise ~2-ulp perturbations of the controls.  Exact arithmetic."""
    key = tuple(pts)
    if key in _kcache:
        return _kcache[key]
    base = O.ExactArc(*pts)
    if base.kind == 'degen':
        _kcache[key] = mp.inf
        return mp.inf
    worst = mp.mpf(0)
    trials = []
    for i in range(3):              # each coordinate alone, both signs
        for j in range(2):
            for sg in (-1, 1):
                q2 = [list(p) for p in pts]
                q2[i][j] = q2[i][j] * (1 + sg * 2.0 ** -51) if q2[i][j] != 0 else sg * 2.0 ** -1074
                trials.append([tuple(p) for p in q2])
    for _ in range(2):
        trials.append([(x * (1 + rng.choice([-1, 1]) * 2.0 ** -51), y * (1 + rng.choice([-1, 1]) * 2.0 ** -51)) for x, y in pts])
    for q2 in trials:
        if base.kind == 'circle':
            q2[2] = q2[0]       # keep CIRCULARSTRING(A,B,A) a full circle under perturbation
        e = O.ExactArc(*q2)
        if e.kind == 'degen' or e.kind != base.kind:
            worst = mp.inf
            break
        worst = max(worst, (mp.sqrt((e.Om[0] - base.Om[0]) ** 2 + (e.Om[1] - base.Om[1]) ** 2) + abs(e.r - base.r)) / base.r)
    _kcache[key] = worst
    return worst


ILL = 1e-9


def check_A(case, res, rng, ora=None):
    a, dev = case
    out = []
    ex = O.ExactArc(*a)
    if res['op'] != 'A':
        return [Finding('exception', 'A.exception', str(res))]
    M = ex.scale()
    if ex.kind == 'degen':
        # no circle exists; proofs decline (valid_arc fails).  Only sanity: GEOS must report linear
        if not res['isLinear'] and not res['isCircle']:
            out.append(Finding('unbacked', 'A.degen_not_linear', 'collinear controls but GEOS isLinear=false'))
        return out
    if res['isLinear']:
        return [Finding('wrong_result', 'A.linear_for_valid_arc',
                        'GEOS treats a valid (non-collinear) arc as a straight line; exact r=%s sweep=%s' % (
                            mp.nstr(ex.r, 12), mp.nstr(ex.sweep(), 12)))]
    # orientation (ArcOrient / Orientation_b64 exact)
    if ex.kind == 'arc' and res['orient'] != ex.orient:
        out.append(Finding('wrong_result', 'A.orientation', 'GEOS orient %d, exact %d' % (res['orient'], ex.orient)))
    # centre / radius
    exact_r = lambda p: (O.ExactArc(*p).r if O.ExactArc(*p).kind != 'degen' else None)
    Sr = perturbed_sensitivity(a, exact_r, rng)
    gr = fnum(res['radius'])
    tol_r = 1e-9 * ex.r + 2.0 ** -40 * M + 64 * Sr
    if not math.isfinite(gr) or abs(gr - ex.r) > tol_r:
        out.append(Finding('wrong_result', 'A.radius', 'GEOS r=%r exact r=%s (tol %s)' % (gr, mp.nstr(ex.r, 17), mp.nstr(tol_r, 3))))
    # sweep / length / area
    th = ex.sweep()
    St = perturbed_sensitivity(a, lambda p: (O.ExactArc(*p).sweep() if O.ExactArc(*p).kind != 'degen' else None), rng)
    ga = fnum(res['angle'])
    tol_t = 1e-9 * th + 1e-12 + 64 * St
    if not math.isfinite(ga) or abs(ga - th) > tol_t:
        cls = 'wrong_result'
        out.append(Finding(cls, 'A.sweep', 'GEOS angle=%r exact sweep=%s (tol %s)' % (ga, mp.nstr(th, 17), mp.nstr(tol_t, 3))))
    L = ex.length()
    gl = fnum(res['length'])
    SL = perturbed_sensitivity(a, lambda p: (O.ExactArc(*p).length() if O.ExactArc(*p).kind != 'degen' else None), rng)
    tol_l = 1e-9 * L + 2.0 ** -40 * M + 64 * SL
    if not math.isfinite(gl) or abs(gl - L) > tol_l:
        out.append(Finding('wrong_result', 'A.length', 'GEOS length=%r exact=%s' % (gl, mp.nstr(L, 17))))
    # proved inequality: chord <= arc length (ArcLength.chord_le_arc_length / ArcChordLength.arc_chord_le_arc_length)
    ch = ex.chord()
    if math.isfinite(gl) and gl < float(ch) * (1 - 1e-12) - 1e-300:
        out.append(Finding('wrong_result', 'A.chord_le_length', 'GEOS length %r < chord %s' % (gl, mp.nstr(ch, 17))))
    S = ex.seg_area()
    gs = fnum(res['area'])
    SA = perturbed_sensitivity(a, lambda p: (O.ExactArc(*p).seg_area() if O.ExactArc(*p).kind != 'degen' else None), rng)
    tol_s = 1e-9 * S + 64 * SA + 1e-300
    if not math.isfinite(gs) or abs(gs - S) > tol_s:
        rel = abs(gs - S) / S if S > 0 else mp.inf
        absrel = abs(gs - S) / (ex.r ** 2 * th + 1e-300)   # error relative to the sector size
        cls = 'wrong_result' if absrel > 1e-9 else 'precision'
        out.append(Finding(cls, 'A.segment_area', 'GEOS area=%r exact=%s rel.err=%s err/(r^2 theta)=%s' % (
            gs, mp.nstr(S, 17), mp.nstr(rel, 3), mp.nstr(absrel, 3))))
    if gs < 0:
        out.append(Finding('wrong_result', 'A.segment_area_nonneg', 'negative segment area %r' % gs))
    # reversal invariance (CurvePolygonOrientation.segment_area_neg magnitude, lin_reverse_commutes)
    rv = res['rev']
    for key in ('radius', 'angle', 'length', 'area'):
        x, y = fnum(res[key]), fnum(rv[key])
        if not (x == y or (math.isnan(x) and math.isnan(y))):
            if abs(x - y) > 1e-12 * max(abs(x), abs(y), 1e-300):
                out.append(Finding('orientation_dependent', 'A.reverse_' + key, '%s fwd %r rev %r' % (key, x, y)))
    # linearization (ArcLinearize / ArcLinearizeBound)
    if 'lin' in res:
        lin = [tuple(map(fnum, p)) for p in res['lin']]
        out += check_lin(ex, a, dev, res, lin, th, M, rng, extra_tol=64 * (Sr + arc_kappa(a, rng) * ex.r))
        rlin = [tuple(map(fnum, p)) for p in rv['lin']]
        if ex.kind == 'circle' or arc_kappa(a, rng) > ILL:
            pass   # CIRCULARSTRING(A,B,A) reversed is the same text; ill-conditioned circles may differ
        elif [a[2]] + lin[:-1] != list(reversed(rlin[:-1] + [a[0]])) and len(lin) == len(rlin):
            fwd = [a[0]] + lin; bwd = list(reversed([a[2]] + rlin))
            dmax = max(math.hypot(p[0] - q2[0], p[1] - q2[1]) for p, q2 in zip(fwd, bwd))
            if dmax > 0:
                out.append(Finding('orientation_dependent', 'A.lin_reverse',
                                   'reversed arc linearizes to different vertices (max diff %.3g)' % dmax))
        elif len(lin) != len(rlin):
            out.append(Finding('orientation_dependent', 'A.lin_reverse_count', 'fwd %d vs rev %d vertices' % (len(lin), len(rlin))))
    # cross-check against the proofs repo's own executable oracle
    if ora is not None and ex.kind == 'arc':
        L_ora, A_ora = ora
        if L_ora not in ('DEGENERATE', 'NAN'):
            lo = parse_hex(L_ora)
            if abs(lo - float(L)) > 1e-12 * float(L) + 2.0 ** -45 * M + 64 * float(SL):
                out.append(Finding('oracle_disagreement', 'ORA.length', 'oracle_arc ARC_LENGTH %r vs mp %s' % (lo, mp.nstr(L, 17))))
            if abs(gl - lo) > tol_l + 1e-12 * lo and not any(f.check == 'A.length' for f in out):
                out.append(Finding('wrong_result', 'A.length_vs_oracle_bin', 'GEOS %r vs oracle_arc %r' % (gl, lo)))
    return out


def check_lin(ex, a, dev, res, lin, th, M, rng, extra_tol=0):
    out = []
    if not lin:
        return [Finding('wrong_result', 'L.empty', 'no linearized points')]
    n_geos = int(res.get('nlin', len(lin)))
    if res.get('linTruncated'):
        return [Finding('wrong_result', 'L.count', 'GEOS produced %d vertices (output truncated); exact sweep %s' % (n_geos, mp.nstr(th, 12)))] \
            if n_geos > 10 * max(2, math.ceil(float(th) / (fnum(res['step']) * math.pi / 180))) else []
    if lin[-1] != a[2]:
        out.append(Finding('wrong_result', 'L.endpoint', 'last vertex %r != p2 %r (ArcLinearize: p_n = gamma(1) exactly)' % (lin[-1], a[2])))
    # vertex count (ArcLinearize.subdiv_n / subdiv_n_ok): n = max(2, ceil(|sweep|/step))
    step = fnum(res['step']) * math.pi / 180
    x = float(th) / step
    n_exp = max(2, math.ceil(x))
    if n_geos != n_exp and abs(x - round(x)) > 1e-6 and not (x < 2 and n_geos == 2):
        out.append(Finding('wrong_result', 'L.count', 'GEOS %d segments, expected max(2,ceil(%.17g/%.17g))=%d' % (n_geos, float(th), step, n_exp)))
    # float evaluation off the exact centre: rounding ~1e-16 relative, far below the 1e-9 tolerance
    ox, oy, r = float(ex.Om[0]), float(ex.Om[1]), float(ex.r)
    pts = [a[0]] + lin
    tol = 1e-9 * r + 2.0 ** -40 * M + float(min(extra_tol, mp.mpf(1e300)))
    stride = max(1, len(pts) // 4000)
    # interior vertices on the exact circle (ArcLinearize.lin_vertex_on_arc)
    worst = 0.0
    for p in lin[:-1:stride]:
        worst = max(worst, abs(math.hypot(p[0] - ox, p[1] - oy) - r))
    if worst > max(tol, 1e-6 * dev):
        out.append(Finding('wrong_result', 'L.on_circle', 'vertex off the exact circle by %.3g (r=%.6g)' % (worst, r)))
    # chord deviation (ArcLinearizeBound.chord_approx_error_bound + chord_sagitta_le_dev)
    maxdev = 0.0
    for i in range(0, len(pts) - 1, stride):
        p, q2 = pts[i], pts[i + 1]
        maxdev = max(maxdev, r - math.hypot((p[0] + q2[0]) / 2 - ox, (p[1] + q2[1]) / 2 - oy))
    if maxdev > dev * (1 + 1e-9) + 2 * tol:
        out.append(Finding('wrong_result', 'L.deviation', 'max chord sagitta %.6g > maxDeviation %r' % (maxdev, dev)))
    # vertices inside the arc span (lin_theta_mono_*, lin_vertex_on_arc)
    if ex.kind == 'arc' and worst <= max(tol, 1e-6 * dev):
        ax, ay, cx, cy = a[0][0], a[0][1], a[2][0], a[2][1]
        Lc = math.hypot(cx - ax, cy - ay)
        bad = 0; eg = None
        for p in lin[:-1:stride]:
            sd = (cx - ax) * (p[1] - ay) - (cy - ay) * (p[0] - ax)
            if abs(sd) > 1e-9 * Lc * (r + M) and (sd > 0) != (ex.sideB > 0):
                bad += 1; eg = eg or p
        if bad:
            out.append(Finding('wrong_result', 'L.span', '%d vertices outside the arc span, e.g. %r' % (bad, eg)))
    if ex.kind == 'circle' and n_geos < 3:
        out.append(Finding('wrong_result', 'L.circle_collapsed', 'full circle linearized to %d vertices' % n_geos))
    return out


def same_pt(p, q2, tol):
    return mpd(p, q2) <= tol


def match_sets(gpts, epts, tol):
    gl = list(gpts); unmatched_e = []
    for e in epts:
        j = next((i for i, g in enumerate(gl) if same_pt(g, e, tol)), None)
        if j is None:
            unmatched_e.append(e)
        else:
            gl.pop(j)
    return unmatched_e, gl


def well_conditioned(info, thr=1e-6):
    if info.get('kappa', 0) > ILL:
        return False
    if info.get('gap_M', 1) < 2.0 ** -36:   # tangency not resolvable in absolute-coordinate doubles
        return False
    return info.get('h_rel', 1) > thr and info.get('span_margin', 1) > thr and info.get('margin', 1) > thr and \
        (info.get('disc_rel', 1) == 0 or abs(info.get('disc_rel', 1)) > thr * thr)


def check_AA(case, res, rng, ora=None):
    a, b = case
    if res['op'] != 'AA':
        return [Finding('exception' if res['op'] == 'EXC' else 'crash', 'AA.exception', str(res))]
    A, B = O.ExactArc(*a), O.ExactArc(*b)
    out = []
    if A.kind == 'degen' or B.kind == 'degen':
        return out
    ab, ba = res['ab'], res['ba']
    for tag, r in (('ab', ab), ('ba', ba)):
        if r['type'] in (1, 2) and r['nPt'] != r['type']:
            out.append(Finding('wrong_result', 'AA.count_type_mismatch', '%s: result type %d but %d points' % (tag, r['type'], r['nPt'])))
    status, E, info = O.arc_arc_points(A, B)
    M = max(A.scale(), B.scale())
    kA, kB = arc_kappa(a, rng), arc_kappa(b, rng)
    info['kappa'] = float(max(kA, kB))
    tol = 1e-9 * max(A.r, B.r) + 2.0 ** -40 * M + 64 * float(min(max(kA * A.r, kB * B.r), mp.mpf(1e300)))
    tol += 2.0 ** -44 * M / max(info.get('h_rel', 1), 1e-300) if 'h_rel' in info else 0
    G = [tuple(map(fnum, p)) for p in ab['pts']]
    if any(not all(map(math.isfinite, p)) for p in G):
        out.append(Finding('wrong_result', 'AA.nan_point', 'non-finite intersection point %r' % G))
    if status == 'coincident':
        # same circle: arcs overlap / touch.  (sanity only: the coincident branch is not proved)
        if max(arc_kappa(a, rng), arc_kappa(b, rng)) > ILL:
            return out
        overl = cocircular_overlap(A, B)
        for tag, r in (('ab', ab), ('ba', ba)):
            for arc in r.get('arcs', []):
                # every reported overlap arc must lie on both inputs: test its three points
                for P in arc:
                    P = tuple(map(fnum, P))
                    dA = O.point_on_arc_dist(A, P); dB = O.point_on_arc_dist(B, P)
                    if max(dA, dB) > 1e-9 * A.r + 2.0 ** -40 * max(A.scale(), B.scale()):
                        out.append(Finding('wrong_result', 'AA.cocircular_arc_off',
                                           '%s: overlap arc %r has point %r off the inputs by %s / %s' % (tag, arc, P, mp.nstr(dA, 4), mp.nstr(dB, 4))))
                        break
        if overl and ab['type'] == 0:
            out.append(Finding('missing_intersection', 'AA.cocircular_missed', 'arcs on one circle overlap/touch (%s) but GEOS reports none' % overl))
        if not overl and ab['type'] != 0:
            out.append(Finding('extra_intersection', 'AA.cocircular_extra', 'disjoint cocircular arcs reported intersecting'))
        return out
    if status == 'concentric':
        if ab['type'] != 0:
            out.append(Finding('extra_intersection', 'AA.concentric', 'concentric distinct circles reported intersecting: %r' % ab))
        return out
    # GEOS may answer COCIRCULAR (type 3) when its inexact centres coincide
    if ab['type'] == 3:
        if well_conditioned(info):
            out.append(Finding('wrong_result', 'AA.false_cocircular', 'distinct circles (exact d^2>0) reported cocircular'))
        return out
    ue, ug = match_sets(G, E, tol + 1e-6 * max(A.r, B.r) * (0 if well_conditioned(info) else 1))
    wc = well_conditioned(info)
    if ue:
        out.append(Finding('missing_intersection' if wc else 'ill_conditioned', 'AA.missing',
                           'expected %s not reported by GEOS (GEOS %r) info=%s' % ([(mp.nstr(x, 17), mp.nstr(y, 17)) for x, y in ue], G, info)))
    if ug:
        # is the extra point actually on both arcs (within tol)?
        for g in ug:
            d1 = O.point_on_arc_dist(A, g); d2 = O.point_on_arc_dist(B, g)
            if max(d1, d2) > tol * 1e3:
                out.append(Finding('extra_intersection' if wc else 'ill_conditioned', 'AA.extra',
                                   'GEOS point %r is %s / %s off the arcs; expected %d points info=%s' % (
                                       g, mp.nstr(d1, 3), mp.nstr(d2, 3), len(E), info)))
    # symmetry (ArcIntersect.arc_arc_intersects_sym)
    Gb = [tuple(map(fnum, p)) for p in ba['pts']]
    u1, u2 = match_sets(Gb, G, tol)
    if (u1 or u2) and wc:
        out.append(Finding('orientation_dependent', 'AA.symmetry', 'intersects(A,B)=%r but intersects(B,A)=%r' % (G, Gb)))
    # distance: 0 iff intersecting (ArcIntersect definitions); closest points on their arcs
    dist = fnum(res['dist'])
    cp = [tuple(map(fnum, p)) for p in res['cp']]
    if wc:
        if E and dist > tol * 10:
            out.append(Finding('wrong_result', 'AA.distance_nonzero', 'arcs intersect at %d points but distance()=%r' % (len(E), dist)))
        if not E:
            dA = O.point_on_arc_dist(A, cp[0]); dB = O.point_on_arc_dist(B, cp[1])
            if max(dA, dB) > tol * 1e3:
                out.append(Finding('wrong_result', 'AA.closest_point_off_arc',
                                   'distance()=%r; closestPoints %r not on the arcs (off by %s, %s)' % (dist, cp, mp.nstr(dA, 4), mp.nstr(dB, 4))))
            elif dist < tol and cp[0] == cp[1]:
                out.append(Finding('wrong_result', 'AA.distance_zero_disjoint', 'disjoint arcs but distance()=%r' % dist))
            # definitional upper bound: endpoints are arc points (unproved composition, sanity only)
            ub = min(math.dist(p, q) for p in (a[0], a[2]) for q in (b[0], b[2]))
            if dist > ub + tol * 10:
                out.append(Finding('wrong_result', 'AA.distance_above_endpoints', 'distance %r exceeds endpoint-pair distance %r' % (dist, ub)))
            # unconditional lower bound (ArcArcDistance.two_circles_dist_lower)
            dO = mp.sqrt((A.Om[0] - B.Om[0]) ** 2 + (A.Om[1] - B.Om[1]) ** 2)
            lb = dO - A.r - B.r
            if dist < lb - tol * 10:
                out.append(Finding('wrong_result', 'AA.distance_lower_bound', 'distance %r below proved bound %s' % (dist, mp.nstr(lb, 12))))
    if ora is not None:
        out += cross_oracle_points('ARC_ARC_XY', ora, E, tol, info)
    return out


def cross_oracle_points(name, ora, E, tol, info):
    toks = ora.split()
    if not toks or toks[0] in ('DEGENERATE', 'NAN', 'COINCIDENT'):
        return []
    n = int(toks[0]); P = [(parse_hex(toks[1 + 2 * i]), parse_hex(toks[2 + 2 * i])) for i in range(n)]
    ue, ug = match_sets(P, E, tol * 10)
    if (ue or ug) and well_conditioned(info):
        return [Finding('oracle_disagreement', 'ORA.' + name, 'oracle_arc %r vs mp %s' % (P, [(mp.nstr(x, 12), mp.nstr(y, 12)) for x, y in E]))]
    return []


def cocircular_overlap(A, B):
    """Do two arcs on the same exact circle share at least one point?  sanity only."""
    if A.kind == 'circle' or B.kind == 'circle':
        return 'full circle'
    for P in (B.A, B.B, B.C):
        if A.contains_exact(P):
            return 'B control in A'
    for P in (A.A, A.B, A.C):
        if B.contains_exact(P):
            return 'A control in B'
    return None


def check_AS(case, res, rng, ora=None):
    a, s = case
    if res['op'] != 'AS':
        return [Finding('exception' if res['op'] == 'EXC' else 'crash', 'AS.exception', str(res))]
    A = O.ExactArc(*a)
    out = []
    if A.kind == 'degen' or s[0] == s[1]:
        return out
    E, info = O.arc_segment_points(A, s[0], s[1])
    M = max(A.scale(), max(abs(c) for p in s for c in p))
    kA = arc_kappa(a, rng); info['kappa'] = float(kA)
    tol = 1e-9 * A.r + 2.0 ** -40 * M + 64 * float(min(kA * A.r, mp.mpf(1e300)))
    tol += 2.0 ** -44 * M / max(info.get('h_rel', 1), 1e-300)
    wc = well_conditioned(info)
    for tag in ('fwd', 'rev'):
        r = res[tag]
        if r['type'] in (1, 2) and r['nPt'] != r['type']:
            out.append(Finding('wrong_result', 'AS.count_type_mismatch', '%s: type %d, %d points' % (tag, r['type'], r['nPt'])))
    G = [tuple(map(fnum, p)) for p in res['fwd']['pts']]
    if any(not all(map(math.isfinite, p)) for p in G):
        out.append(Finding('wrong_result', 'AS.nan_point', 'non-finite point %r' % G))
    ue, ug = match_sets(G, E, tol)
    if ue:
        exact_end = any(isinstance(e[0], mp.mpf) and any(abs(e[0] - p[0]) == 0 and abs(e[1] - p[1]) == 0 for p in s) for e in ue)
        out.append(Finding('missing_intersection' if (wc or (exact_end and kA <= ILL)) else 'ill_conditioned', 'AS.missing',
                           'expected %s not reported (GEOS %r) info=%s%s' % (
                               [(mp.nstr(x, 17), mp.nstr(y, 17)) for x, y in ue], G, info, ' [segment endpoint exactly on arc]' if exact_end else '')))
    for g in ug:
        dA = O.point_on_arc_dist(A, g); dS = O.seg_point_distance(s[0], s[1], g)
        if max(dA, dS) > tol * 1e3:
            out.append(Finding('extra_intersection' if wc else 'ill_conditioned', 'AS.extra',
                               'GEOS point %r is %s off arc, %s off segment' % (g, mp.nstr(dA, 3), mp.nstr(dS, 3))))
        elif max(dA, dS) > tol * 10 and wc:
            out.append(Finding('precision', 'AS.inaccurate_point', 'GEOS point %r is %s off arc, %s off segment (tol %s)' % (
                g, mp.nstr(dA, 3), mp.nstr(dS, 3), mp.nstr(tol, 3))))
    Gr = [tuple(map(fnum, p)) for p in res['rev']['pts']]
    u1, u2 = match_sets(Gr, G, tol)
    if (u1 or u2) and wc:
        out.append(Finding('orientation_dependent', 'AS.segment_direction', 'segment P->Q gives %r, Q->P gives %r' % (G, Gr)))
    dist = fnum(res['dist'])
    cp = [tuple(map(fnum, p)) for p in res['cp']]
    if wc:
        if E and dist > tol * 10:
            out.append(Finding('wrong_result', 'AS.distance_nonzero', 'intersecting but distance()=%r' % dist))
        if not E:
            dA = O.point_on_arc_dist(A, cp[0]); dS = O.seg_point_distance(s[0], s[1], cp[1])
            dA2 = O.point_on_arc_dist(A, cp[1]); dS2 = O.seg_point_distance(s[0], s[1], cp[0])
            if max(dA, dS) > tol * 1e3 and max(dA2, dS2) <= tol * 1e3:
                out.append(Finding('wrong_result', 'AS.closest_points_swapped', 'closestPoints returned (segment point, arc point) %r; distance()=%r' % (cp, dist)))
            elif max(dA, dS) > tol * 1e3:
                out.append(Finding('wrong_result', 'AS.closest_point_off', 'closestPoints %r off arc by %s / off segment by %s; distance()=%r' % (
                    cp, mp.nstr(dA, 4), mp.nstr(dS, 4), dist)))
    if ora is not None:
        out += cross_oracle_points('ARC_SEGMENT_XY', ora, E, tol, info)
    return out


def check_AP(case, res, rng, ora=None):
    a, p = case
    if res['op'] != 'AP':
        return [Finding('exception' if res['op'] == 'EXC' else 'crash', 'AP.exception', str(res))]
    A = O.ExactArc(*a)
    if A.kind == 'degen':
        return []
    out = []
    dE = O.arc_point_distance(A, p)        # LECArcRow.arc_dist_exact
    M = max(A.scale(), abs(p[0]), abs(p[1]))
    S = perturbed_sensitivity(list(a) + [p], lambda q: (O.arc_point_distance(O.ExactArc(*q[:3]), q[3]) if O.ExactArc(*q[:3]).kind != 'degen' else None), rng)
    tol = 1e-9 * max(dE, A.r) + 2.0 ** -40 * M + 64 * S
    g = fnum(res['dist'])
    if not math.isfinite(g) or abs(g - dE) > tol:
        out.append(Finding('wrong_result', 'AP.distance', 'GEOS %r exact %s (tol %s)' % (g, mp.nstr(dE, 17), mp.nstr(tol, 3))))
    if A.kind == 'arc' or A.kind == 'circle':
        P = O.q(p)
        on = A.contains_exact(P)
        if on != res['contains'] and arc_kappa(a, rng) <= ILL:
            out.append(Finding('wrong_result', 'AP.containsPoint' + ('_fullcircle' if A.kind == 'circle' else ''), 'GEOS containsPoint=%s, exact %s' % (res['contains'], on)))
    if ora is not None and ora not in ('DEGENERATE', 'NAN'):
        o = parse_hex(ora)
        if abs(o - float(dE)) > tol * 10:
            out.append(Finding('oracle_disagreement', 'ORA.ARC_DISTANCE', 'oracle_arc %r vs exact %s' % (o, mp.nstr(dE, 17))))
    return out

# ---- geometry-level checks ----

def parse_wkt_curve(w):
    """Tiny parser for the WKT we generate: returns list of (kind, pts) sections and polygon flag."""
    import re
    poly = w.startswith('CURVEPOLYGON')
    secs = []
    for m in re.finditer(r'(CIRCULARSTRING|LINESTRING)?\s*\(([-0-9eE+.,\s]+)\)', w):
        kind = 'arc' if m.group(1) == 'CIRCULARSTRING' else 'line'
        pts = [tuple(map(float, t.split())) for t in m.group(2).split(',')]
        secs.append((kind, pts))
    return poly, secs


def sections_exact(secs):
    L = mp.mpf(0); A2 = mp.mpf(0); ok = True
    for kind, pts in secs:
        if kind == 'arc':
            for i in range(0, len(pts) - 2, 2):
                ex = O.ExactArc(pts[i], pts[i + 1], pts[i + 2])
                if ex.kind == 'degen':
                    L += mpd(pts[i], (mp.mpf(pts[i + 2][0]), mp.mpf(pts[i + 2][1])))
                    ok = False
                    continue
                L += ex.length()
                p0, p2 = pts[i], pts[i + 2]
                A2 += (mp.mpf(p0[0]) * p2[1] - mp.mpf(p2[0]) * p0[1])
                sgnv = 1 if ex.orient > 0 else -1     # GEOS full circle: COLLINEAR -> treated as CW
                A2 += 2 * sgnv * ex.seg_area()
        else:
            for p, q2 in zip(pts, pts[1:]):
                L += mpd(p, (mp.mpf(q2[0]), mp.mpf(q2[1])))
                A2 += (mp.mpf(p[0]) * q2[1] - mp.mpf(q2[0]) * p[1])
    return L, abs(A2) / 2, ok


def check_G(line, res, rng, ora=None):
    if res['op'] != 'G':
        return [Finding('exception' if res['op'] == 'EXC' else 'crash', 'G.exception', str(res))]
    w = line[2:]
    poly, secs = parse_wkt_curve(w)
    L, A, ok = sections_exact(secs)
    out = []
    M = max(abs(c) for _, pts in secs for p in pts for c in p) + 1e-300
    gl, ga = fnum(res['length']), fnum(res['area'])
    if not ok:
        return out
    kmax = max([arc_kappa(tuple(pts[i:i + 3]), rng) for k, pts in secs if k == 'arc' for i in range(0, len(pts) - 2, 2)] + [mp.mpf(0)])
    if kmax > ILL:
        return [Finding('ill_conditioned', 'G.ill', 'an arc has circle conditioning %s' % mp.nstr(kmax, 3))] if False else []
    tol = 1e-9 * L + 2.0 ** -40 * M * len(secs[0][1])
    if not math.isfinite(gl) or abs(gl - L) > tol:
        out.append(Finding('wrong_result', 'G.length', 'GEOS length %r, exact sum of arc lengths %s (CurveLength.curve_length_additive)' % (gl, mp.nstr(L, 17))))
    if poly:
        tA = 1e-9 * A + 2.0 ** -40 * M * M
        if not math.isfinite(ga) or abs(ga - A) > tA:
            out.append(Finding('wrong_result', 'G.area', 'GEOS area %r, shoelace+segments %s (composition unproved: sanity)' % (ga, mp.nstr(A, 17))))
        if res['linType'] != 'Polygon':
            out.append(Finding('wrong_type', 'G.lin_type', 'linearized CurvePolygon is %s' % res['linType']))
    else:
        if res['linType'] != 'LineString':
            out.append(Finding('wrong_type', 'G.lin_type', 'linearized curve is %s' % res['linType']))
    # inscribed polyline is no longer than the curve (chord_le_arc_length per piece)
    ll = fnum(res['linLength'])
    if ll > L * (1 + 1e-9) + tol:
        out.append(Finding('wrong_result', 'G.lin_longer', 'linearized length %r > curve length %s' % (ll, mp.nstr(L, 17))))
    return out


def out_piece_length(piece):
    t = piece['t']
    if t == 'CompoundCurve':
        return sum(out_piece_length(p) for p in piece['parts']), \
            sum(out_piece_length(p)[1] for p in piece['parts']) if False else None
    return None


def piece_lengths(piece):
    """(curved length, straight length) of one noded output piece, exact on its own coords."""
    t = piece['t']
    if t == 'CompoundCurve':
        c = s = mp.mpf(0)
        for p in piece['parts']:
            a, b = piece_lengths(p); c += a; s += b
        return c, s
    pts = [tuple(map(fnum, p)) for p in piece['c']]
    if t == 'CircularString':
        L = mp.mpf(0)
        for i in range(0, len(pts) - 2, 2):
            ex = O.ExactArc(pts[i], pts[i + 1], pts[i + 2])
            L += ex.length() if ex.kind != 'degen' else mpd(pts[i], (mp.mpf(pts[i + 2][0]), mp.mpf(pts[i + 2][1])))
        return L, mp.mpf(0)
    return mp.mpf(0), sum((mpd(p, (mp.mpf(q2[0]), mp.mpf(q2[1]))) for p, q2 in zip(pts, pts[1:])), mp.mpf(0))


def piece_endpoints(piece):
    if piece['t'] == 'CompoundCurve':
        e = []
        for p in piece['parts']:
            e += piece_endpoints(p)
        return e
    c = piece['c']
    return [tuple(map(fnum, c[0])), tuple(map(fnum, c[-1]))]


def check_N(line, res, rng, ora=None):
    if res['op'] != 'N':
        return [Finding('exception' if res['op'] == 'EXC' else 'crash', 'N.exception', str(res))]
    w = line.split(None, 2)[2]
    inner = w[len('GEOMETRYCOLLECTION('):-1]
    # split the two members at the top-level comma
    depth = 0; cut = None
    for i, ch in enumerate(inner):
        depth += ch == '('; depth -= ch == ')'
        if ch == ',' and depth == 0:
            cut = i; break
    g1, g2 = inner[:cut].strip(), inner[cut + 1:].strip()
    _, s1 = parse_wkt_curve(g1); _, s2 = parse_wkt_curve(g2)
    out = []
    # input lengths split into curved and straight parts
    def lens(secs):
        c = s = mp.mpf(0); ok = True
        for kind, pts in secs:
            if kind == 'arc':
                for i in range(0, len(pts) - 2, 2):
                    ex = O.ExactArc(pts[i], pts[i + 1], pts[i + 2])
                    if ex.kind == 'degen':
                        ok = False
                    else:
                        c += ex.length()
            else:
                for p, q2 in zip(pts, pts[1:]):
                    s += mpd(p, (mp.mpf(q2[0]), mp.mpf(q2[1])))
        return c, s, ok
    c1, l1, ok1 = lens(s1); c2, l2, ok2 = lens(s2)
    if not (ok1 and ok2):
        return out
    arcs_all = [tuple(pts[i:i + 3]) for secs in (s1, s2) for k, pts in secs if k == 'arc' for i in range(0, len(pts) - 2, 2)]
    if any(arc_kappa(x, rng) > ILL for x in arcs_all):
        return out
    overlap = False
    A1 = [O.ExactArc(*x) for k, pts in s1 if k == 'arc' for x in [tuple(pts[i:i + 3]) for i in range(0, len(pts) - 2, 2)]]
    A2 = [O.ExactArc(*x) for k, pts in s2 if k == 'arc' for x in [tuple(pts[i:i + 3]) for i in range(0, len(pts) - 2, 2)]]
    for x in A1:
        for y in A2:
            if x.O == y.O and x.r2 == y.r2:
                overlap = True
    segs = lambda secs: set(frozenset((p, q2)) for k, pts in secs if k == 'line' for p, q2 in zip(pts, pts[1:]))
    if segs(s1) & segs(s2):
        overlap = True       # identical straight edges are deduplicated by design
    oc = os_ = mp.mpf(0)
    ends = []
    slivers = []
    ends_curved, ends_straight = [], []
    for piece in res['out']:
        a, b = piece_lengths(piece); oc += a; os_ += b
        ends += piece_endpoints(piece)
        for part in (piece['parts'] if piece['t'] == 'CompoundCurve' else [piece]):
            (ends_curved if part['t'] == 'CircularString' else ends_straight).extend(piece_endpoints(part))
        for part in (piece['parts'] if piece['t'] == 'CompoundCurve' else [piece]):
            pa, pb = piece_lengths(part)
            if 0 < pa + pb < 1e-12 * (c1 + c2 + l1 + l2):
                slivers.append((part['c'], pa + pb))
    M = max(abs(c) for secs in (s1, s2) for _, pts in secs for p in pts for c in p) + 1e-300
    Rmax = max([O.ExactArc(*pts[i:i + 3]).r for secs in (s1, s2) for k, pts in secs if k == 'arc'
                for i in range(0, len(pts) - 2, 2)] + [mp.mpf(1e-300)])
    tol = 1e-7 * (c1 + c2 + l1 + l2) + 2.0 ** -36 * M
    # noding partitions the input point set (ArcSplitAtNode.arc_split_cover / arc_split_disjoint,
    # CurveLength.curve_length_additive): total length is preserved
    if not overlap and abs((oc + os_) - (c1 + c2 + l1 + l2)) > tol:
        out.append(Finding('wrong_result', 'N.length_not_preserved', 'input length %s, noded output length %s' % (
            mp.nstr(c1 + c2 + l1 + l2, 15), mp.nstr(oc + os_, 15))))
    # pieces far below the input scale are rounding slivers (sanity only: not a proved property)
    if slivers and not any(min_input_piece < 1e-11 * (c1 + c2 + l1 + l2) for min_input_piece in [min_piece_len(s1 + s2)]):
        out.append(Finding('wrong_result', 'N.sliver_piece', '%d output piece(s) of length ~%s, e.g. %r' % (
            len(slivers), mp.nstr(slivers[0][1], 3), slivers[0][0])))
    # arcs must stay arcs (the curve-preserving noder is the point of GeometryNoder on curves)
    if (c1 + c2) > 0 and os_ > (l1 + l2) + tol:
        out.append(Finding('wrong_type', 'N.linearized', 'straight input length %s but %s straight in output (arcs became straight)' % (
            mp.nstr(l1 + l2, 12), mp.nstr(os_, 12))))
    # every well-conditioned crossing must be a node
    exp_nodes = []
    for k1, p1 in s1:
        for k2, p2 in s2:
            exp_nodes += pair_nodes(k1, p1, k2, p2)
    tolp = 1e-9 * Rmax + 2.0 ** -40 * M
    for X, wc, aends, sends in exp_nodes:
        if not wc or overlap:
            continue
        near = lambda pts: any(mpd(e, X) <= tolp * 1e3 for e in pts)
        if not near(ends):
            out.append(Finding('missing_intersection', 'N.missing_node', 'crossing (%s, %s) is not a node of the noded output' % (
                mp.nstr(X[0], 17), mp.nstr(X[1], 17))))
        elif not near(aends) and not near(ends_curved):
            out.append(Finding('missing_intersection', 'N.arc_not_split', 'node (%s, %s) lies inside an input arc but no curved piece ends there' % (
                mp.nstr(X[0], 17), mp.nstr(X[1], 17))))
        elif sends is not None and not near(sends) and not near(ends_straight):
            out.append(Finding('missing_intersection', 'N.line_not_split', 'node (%s, %s) lies inside an input segment but no straight piece ends there' % (
                mp.nstr(X[0], 17), mp.nstr(X[1], 17))))
    return out


def min_piece_len(secs):
    m = mp.inf
    for k, pts in secs:
        for p, q2 in zip(pts, pts[1:]):
            m = min(m, mpd(p, (mp.mpf(q2[0]), mp.mpf(q2[1]))))
    return m


def pair_nodes(k1, p1, k2, p2):
    """Expected nodes between two input sections: (X, well_conditioned, arcs_ends, seg_ends)
    where arcs_ends / seg_ends are the input endpoints X may coincide with (if X is
    interior to an arc it must become an endpoint of a curved piece, likewise for lines)."""
    res = []
    arcs1 = [p1[i:i + 3] for i in range(0, len(p1) - 2, 2)] if k1 == 'arc' else []
    arcs2 = [p2[i:i + 3] for i in range(0, len(p2) - 2, 2)] if k2 == 'arc' else []
    segs1 = list(zip(p1, p1[1:])) if k1 == 'line' else []
    segs2 = list(zip(p2, p2[1:])) if k2 == 'line' else []
    for a in arcs1:
        A = O.ExactArc(*a)
        if A.kind == 'degen':
            continue
        for b in arcs2:
            B = O.ExactArc(*b)
            if B.kind == 'degen':
                continue
            st, E, info = O.arc_arc_points(A, B)
            if st == 'ok':
                res += [(X, well_conditioned(info), [a[0], a[2], b[0], b[2]], None) for X in E]
        for sg in segs2:
            E, info = O.arc_segment_points(A, sg[0], sg[1])
            res += [(X, well_conditioned(info), [a[0], a[2]], [sg[0], sg[1]]) for X in E]
    for sg in segs1:
        for b in arcs2:
            B = O.ExactArc(*b)
            if B.kind == 'degen':
                continue
            E, info = O.arc_segment_points(B, sg[0], sg[1])
            res += [(X, well_conditioned(info), [b[0], b[2]], [sg[0], sg[1]]) for X in E]
    return res


def check_SS(case, res, rng, ora=None):
    r = res['r']
    if r['type'] in (1, 2) and r['nPt'] != r['type']:
        return [Finding('wrong_result', 'SS.count_type_mismatch', 'result type %d but getNumPoints()=%d' % (r['type'], r['nPt']))]
    return []


CHECKS = {'SS': check_SS, 'A': check_A, 'AA': check_AA, 'AS': check_AS, 'AP': check_AP, 'G': check_G, 'N': check_N}


def oracle_queries(op, data):
    if op == 'A':
        a, _ = data
        return [('ARC_LENGTH', a), ('ARC_AREA', a)]
    if op == 'AA':
        a, b = data
        return [('ARC_ARC_XY', a + b)]
    if op == 'AS':
        a, s = data
        return [('ARC_SEGMENT_XY', a + s)]
    if op == 'AP':
        a, p = data
        return [('ARC_DISTANCE', a + [p])]
    return []


def evaluate(lines, rng, use_oracle=True):
    results = run_probe(lines)
    parsed = [parse_case(l) for l in lines]
    ora_out = {}
    if use_oracle:
        qs = []; owner = []
        for i, (op, data) in enumerate(parsed):
            for qq in oracle_queries(op, data):
                qs.append(qq); owner.append(i)
        if qs:
            outs = run_oracle(qs)
            for i, o in zip(owner, outs):
                ora_out.setdefault(i, []).append(o)
    findings = []
    for i, ((op, data), res) in enumerate(zip(parsed, results)):
        if res.get('op') == 'CRASH':
            findings.append([Finding('crash', op + '.crash', str(res))]); continue
        o = ora_out.get(i)
        try:
            if op == 'A':
                f = check_A(data, res, rng, (o[0], o[1]) if o else None)
            elif op in ('AA', 'AS', 'AP'):
                f = CHECKS[op](data, res, rng, o[0] if o else None)
            else:
                f = CHECKS[op](data, res, rng)
        except Exception as e:  # harness bug: report, do not hide
            f = [Finding('harness_error', op + '.harness', repr(e))]
        findings.append(f)
    return results, findings

# ---------------------------------------------------------------------------
# shrinking
# ---------------------------------------------------------------------------

def signature(fs, keep=None):
    s = sorted(set(f.check for f in fs if f.cls not in ('ill_conditioned', 'oracle_disagreement', 'unbacked')))
    return tuple(s if keep is None else [c for c in s if c in keep])


def fails(line, target, rng):
    try:
        _, f = evaluate([line], rng, use_oracle=False)
    except Exception:
        return False
    return set(target) <= set(signature(f[0]))


def shrink(line, rng, budget=400):
    _, f = evaluate([line], rng, use_oracle=False)
    target = signature(f[0])
    if not target:
        return line, target
    target = target[:1]
    op = line.split()[0]
    if op in ('G', 'N'):
        return line, target
    head = op
    nums = line.split()[1:]
    tries = 0
    best = nums[:]
    improved = True
    while improved and tries < budget:
        improved = False
        for i in range(len(best)):
            v = float(best[i])
            cands = []
            for k in range(0, 17):
                c = round(v, k)
                if c != v:
                    cands.append(c)
            cands = [0.0] + [float(int(v))] + cands
            seen = set()
            for c in cands:
                tries += 1
                rc = repr(float(c)).replace('.0', '') if float(c).is_integer() and abs(c) < 1e15 else repr(float(c))
                if rc in seen or rc == best[i]:
                    continue
                seen.add(rc)
                trial = best[:]; trial[i] = rc
                if fails(head + ' ' + ' '.join(trial), target, rng):
                    if len(rc) < len(best[i]) or float(rc) == 0:
                        best = trial; improved = True
                        break
    return head + ' ' + ' '.join(best), target

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def cmd_run(args):
    rng = random.Random(args.seed)
    chk_rng = random.Random(args.seed + 1)
    gen = Gen(rng)
    t_end = time.time() + 60 * args.minutes
    os.makedirs(args.out, exist_ok=True)
    counts = collections.Counter(); fcount = collections.Counter(); cls_count = collections.Counter()
    examples = collections.defaultdict(list)
    n = 0; batch = args.batch
    log = open(os.path.join(args.out, 'findings.jsonl'), 'w')
    while time.time() < t_end and n < args.max_cases:
        lines, tags = [], []
        for _ in range(batch):
            l, t = gen.case(); lines.append(l); tags.append(t)
        _, fs = evaluate(lines, chk_rng, use_oracle=not args.no_oracle)
        for idx, (l, t, f) in enumerate(zip(lines, tags, fs)):
            counts[t] += 1
            for x in f:
                fcount[(x.cls, x.check)] += 1; cls_count[x.cls] += 1
                if len(examples[(x.cls, x.check)]) < 5:
                    examples[(x.cls, x.check)].append({'case': l, 'gen': t, 'msg': x.msg})
                log.write(json.dumps({'i': n + idx, 'gen': t, 'case': l, 'cls': x.cls, 'check': x.check, 'msg': x.msg}) + '\n')
        n += len(lines)
        if args.verbose:
            print('%d cases, %d finding-instances, %.0fs left' % (n, sum(fcount.values()), t_end - time.time()), flush=True)
    log.close()
    summary = {'seed': args.seed, 'cases': n, 'by_generator': dict(counts),
               'findings': {'%s|%s' % k: v for k, v in sorted(fcount.items())},
               'by_class': dict(cls_count),
               'examples': {'%s|%s' % k: v for k, v in examples.items()}}
    json.dump(summary, open(os.path.join(args.out, 'summary.json'), 'w'), indent=1)
    print(json.dumps({k: summary[k] for k in ('seed', 'cases', 'by_class', 'findings')}, indent=1))


def cmd_case(args):
    rng = random.Random(1)
    res, fs = evaluate([args.line], rng)
    print(json.dumps(res[0])[:4000])
    for f in fs[0]:
        print('%-22s %-28s %s' % (f.cls, f.check, f.msg))
    if not fs[0]:
        print('no findings')


def cmd_shrink(args):
    rng = random.Random(1)
    l, t = shrink(args.line, rng)
    print(t)
    print(l)


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--seed', type=int, default=1); r.add_argument('--minutes', type=float, default=40)
    r.add_argument('--max-cases', type=int, default=10 ** 9); r.add_argument('--batch', type=int, default=500)
    r.add_argument('--out', default=os.path.join(HERE, 'results')); r.add_argument('--no-oracle', action='store_true')
    r.add_argument('-v', '--verbose', action='store_true')
    c = sp.add_parser('case'); c.add_argument('line')
    s = sp.add_parser('shrink'); s.add_argument('line')
    a = ap.parse_args()
    {'run': cmd_run, 'case': cmd_case, 'shrink': cmd_shrink}[a.cmd](a)


if __name__ == '__main__':
    main()
