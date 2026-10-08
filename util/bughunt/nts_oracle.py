"""Exact / high-precision oracle for GEOS circular arcs, transcribed from
NetTopologySuite.Proofs (github.com/grootstebozewolf/NetTopologySuite.Proofs,
pinned commit 57f82369) and cross-checked at run time against that repo's own
reference executable (oracle/driver.ml + the Coq-extracted kernel, built by
build_oracle.sh as bin/oracle_arc).

Every function cites the theorem / file it transcribes.  Rationals are exact
(fractions.Fraction: every binary64 is dyadic, as in driver.ml's zarith Q
layer).  Where a square root or trig function is unavoidable we use mpmath at
60 significant digits, applied to the exact rational invariants (the same
"one rounding past the certified algebra" discipline as driver.ml).
"""
from fractions import Fraction as Fr
import math
import mpmath as mp

mp.mp.dps = 60
PI = mp.pi


def q(p):
    return (Fr(p[0]), Fr(p[1]))


def sgn(v):
    return (v > 0) - (v < 0)


def cross(o, a, b):
    """Twice the signed triangle area; theories/ArcOrient.v cross_R_pt,
    theories/CurveGeometry.v valid_arc.  Exact."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def circumcentre(a, b, c):
    """theories/CurveGeometry.v arc_center (Graphics-Gems / determinant form),
    proved equidistant in theories/CircumcentreQSound.v
    circumcentre_formula_equidistant.  Same algebra as driver.ml circumcentre_q."""
    ax, ay = a; bx, by = b; cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if d == 0:
        return None
    na = ax * ax + ay * ay; nb = bx * bx + by * by; nc = cx * cx + cy * cy
    ox = (na * (by - cy) + nb * (cy - ay) + nc * (ay - by)) / d
    oy = (na * (cx - bx) + nb * (ax - cx) + nc * (bx - ax)) / d
    return (ox, oy)


def mpf(x):
    return mp.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else mp.mpf(x)


class ExactArc:
    """One SQL/MM three-point arc.

    kind: 'arc'    valid_arc holds (non-collinear controls)
          'circle' p0 == p2 != p1: GEOS/ADR-0005 full circle, centre = midpoint(p0,p1)
                   (theories/IntakeWalkerCircle.v: CIRCULARSTRING(A,B,A) == CIRCLE(A,B,ogc_c))
          'degen'  collinear controls or p0 == p1 == p2 (no circle; proofs decline)
    """

    def __init__(self, p0, p1, p2):
        self.f = (tuple(map(float, p0)), tuple(map(float, p1)), tuple(map(float, p2)))
        self.A, self.B, self.C = q(p0), q(p1), q(p2)
        A, B, C = self.A, self.B, self.C
        if A == C:
            if A == B:
                self.kind = 'degen'
                return
            self.kind = 'circle'
            self.O = ((A[0] + B[0]) / 2, (A[1] + B[1]) / 2)
            self.r2 = ((B[0] - A[0]) ** 2 + (B[1] - A[1]) ** 2) / 4
            self.orient = 0
            self.major = 1
            self.cos = Fr(1)
        else:
            D = cross(A, B, C)
            if D == 0:
                self.kind = 'degen'
                return
            self.kind = 'arc'
            self.orient = sgn(D)          # +1 CCW, -1 CW (Orientation index convention)
            self.O = circumcentre(A, B, C)
            O = self.O
            vax, vay = A[0] - O[0], A[1] - O[1]
            vcx, vcy = C[0] - O[0], C[1] - O[1]
            self.r2 = vax * vax + vay * vay
            # driver.ml arc_invariants_q: cos Theta0 (AngleBetween.cos_angle_between)
            self.cos = (vax * vcx + vay * vcy) / self.r2
            sb = sgn(self.side(B)); so = sgn(self.side(O))
            self.major = 1 if (so != 0 and sb == so) else 0
            self.sideB = sb
        self.Om = (mpf(self.O[0]), mpf(self.O[1]))
        self.r = mp.sqrt(mpf(self.r2))

    # chord-side sign test (theories/LECArcRow.v arc_span_contains_iff_sign,
    # theories/ArcSpanAtan2.v arc_span_contains_atan2_iff_chord_sign)
    def side(self, X):
        A, C = self.A, self.C
        return (C[0] - A[0]) * (X[1] - A[1]) - (C[1] - A[1]) * (X[0] - A[0])

    def side_mp(self, X):
        ax, ay = mpf(self.A[0]), mpf(self.A[1])
        cx, cy = mpf(self.C[0]), mpf(self.C[1])
        return (cx - ax) * (X[1] - ay) - (cy - ay) * (X[0] - ax)

    def in_span_mp(self, X):
        """X assumed on the circle.  Closed sign test 0 <= side(B)*side(X)."""
        if self.kind == 'circle':
            return True
        s = self.side_mp(X)
        scale = (abs(mpf(self.C[0] - self.A[0])) + abs(mpf(self.C[1] - self.A[1]))) * (self.r + 1)
        if abs(s) <= scale * mp.mpf(10) ** (-45):
            return True       # on the chord line => an endpoint (on_circle_side_zero_is_endpoint)
        return sgn(s) == self.sideB

    def span_margin(self, X):
        """Relative distance of X from the span boundary (for conditioning)."""
        if self.kind == 'circle':
            return mp.inf
        dA = mp.sqrt((X[0] - mpf(self.A[0])) ** 2 + (X[1] - mpf(self.A[1])) ** 2)
        dC = mp.sqrt((X[0] - mpf(self.C[0])) ** 2 + (X[1] - mpf(self.C[1])) ** 2)
        return min(dA, dC) / self.r

    def sweep(self):
        """|Theta|: driver.ml run_arc_length (exact (1-cos)/2 then asin),
        theories/ArcLength.v chord_subtended_sq; full circle 2*pi
        (ArcLength.arc_length_full_turn)."""
        if self.kind == 'circle':
            return 2 * PI
        s = mp.sqrt(mpf((1 - self.cos) / 2))
        t0 = 2 * mp.asin(s)
        return 2 * PI - t0 if self.major else t0

    def length(self):
        """theories/ArcLength.v arc_length r*theta."""
        return self.r * self.sweep()

    def seg_area(self):
        """theories/ArcArea.v segment_area (r^2/2)(theta - sin theta);
        full disc pi r^2 (segment_area_full_disc)."""
        th = self.sweep()
        return mpf(self.r2) / 2 * (th - mp.sin(th))

    def chord(self):
        return mp.sqrt(mpf((self.C[0] - self.A[0]) ** 2 + (self.C[1] - self.A[1]) ** 2))

    def scale(self):
        return max(abs(c) for p in self.f for c in p) + 1e-300

    def on_circle_exact(self, P):
        """inCircle_R == 0 decided exactly (theories/ArcArcCircles.v
        inCircle_R_zero_implies_equidistant)."""
        return (P[0] - self.O[0]) ** 2 + (P[1] - self.O[1]) ** 2 == self.r2

    def contains_exact(self, P):
        """Exact point-on-arc for a rational point (ArcIntersect.arc_span_contains)."""
        if self.kind == 'degen' or not self.on_circle_exact(P):
            return False
        if self.kind == 'circle' or P == self.A or P == self.C:
            return True
        return sgn(self.side(P)) == self.sideB


def angle_of(arc, X):
    return mp.atan2(X[1] - arc.Om[1], X[0] - arc.Om[0])


# ---------------------------------------------------------------------------
# Arc / arc intersection
# ---------------------------------------------------------------------------

def arc_arc_points(a1, a2):
    """Expected intersection point set of two non-degenerate arcs.

    Returns (status, points, info).  status in {'ok','coincident','concentric'}.
    Circle-circle points: radical-line construction, theories/ArcArcCircles.v
    radical_points_on_circles (both points on both circles) and
    two_circles_radical_point_unique (no other common points).  Tangency is
    decided exactly in Q (h^2 = r1^2 - a^2 with a = (d^2+r1^2-r2^2)/(2d), so
    4 d^2 h^2 = 4 d^2 r1^2 - (d^2+r1^2-r2^2)^2 is rational).  Span membership:
    LECArcRow.arc_span_contains_iff_sign.  Shared endpoints are exact
    (ArcArcSound.arc_arc_intersects_shared_vertex)."""
    O1, O2 = a1.O, a2.O
    dx, dy = O2[0] - O1[0], O2[1] - O1[1]
    d2 = dx * dx + dy * dy
    if d2 == 0:
        if a1.r2 == a2.r2:
            return 'coincident', [], {}
        return 'concentric', [], {}
    k = d2 + a1.r2 - a2.r2
    disc4 = 4 * d2 * a1.r2 - k * k          # = 4 d^2 h^2
    info = {'disc_rel': float(mpf(disc4) / (4 * mpf(d2) * mpf(a1.r2)))}
    d = mp.sqrt(mpf(d2))
    M = max(a1.scale(), a2.scale())
    info['gap_M'] = float(min(abs(d - (a1.r + a2.r)), abs(d - abs(a1.r - a2.r))) / M)
    if disc4 < 0:
        return 'ok', [], info
    a = mpf(k) / (2 * d)
    h = mp.sqrt(mpf(disc4) / (4 * mpf(d2)))
    ux, uy = mpf(dx) / d, mpf(dy) / d
    mx, my = a1.Om[0] + a * ux, a1.Om[1] + a * uy
    cands = [(mx - h * uy, my + h * ux)]
    if disc4 != 0:
        cands.append((mx + h * uy, my - h * ux))
    pts = []
    margins = []
    for X in cands:
        # exact snap to a shared control endpoint when it is a true common point
        snapped = None
        for P in (a1.A, a1.C, a2.A, a2.C):
            if a1.on_circle_exact(P) and a2.on_circle_exact(P):
                Pm = (mpf(P[0]), mpf(P[1]))
                if mp.sqrt((Pm[0] - X[0]) ** 2 + (Pm[1] - X[1]) ** 2) <= (a1.r + a2.r) * mp.mpf(10) ** -40:
                    snapped = P
        if snapped is not None:
            if a1.contains_exact(snapped) and a2.contains_exact(snapped):
                pts.append((mpf(snapped[0]), mpf(snapped[1])))
            continue
        if a1.in_span_mp(X) and a2.in_span_mp(X):
            pts.append(X)
        margins.append(min(a1.span_margin(X), a2.span_margin(X)))
    info['h_rel'] = float(h / min(a1.r, a2.r))
    info['span_margin'] = float(min(margins)) if margins else float('inf')
    return 'ok', pts, info


# ---------------------------------------------------------------------------
# Arc / segment intersection
# ---------------------------------------------------------------------------

def arc_segment_points(arc, P, Q):
    """Circle-line points F +- h u (theories/ArcSegmentCircles.v
    line_circle_radical_point, arc_line_circle_intersect), kept iff in the arc
    span (LECArcRow.arc_span_contains_iff_sign) and on the segment
    (ArcIntersect.arc_chord_intersects: between P Q X).  Foot parameter and
    h^2 are exact rationals; segment endpoints lying exactly on the arc are
    detected exactly."""
    P, Q = q(P), q(Q)
    O = arc.O
    dx, dy = Q[0] - P[0], Q[1] - P[1]
    L2 = dx * dx + dy * dy
    if L2 == 0:
        return [], {}
    s = ((O[0] - P[0]) * dx + (O[1] - P[1]) * dy) / L2
    F = (P[0] + s * dx, P[1] + s * dy)
    h2 = arc.r2 - ((O[0] - F[0]) ** 2 + (O[1] - F[1]) ** 2)
    info = {'h_rel': float(mp.sqrt(abs(mpf(h2))) / arc.r)}
    dl = mp.sqrt((mpf(O[0]) - mpf(F[0])) ** 2 + (mpf(O[1]) - mpf(F[1])) ** 2)
    M = max(arc.scale(), max(abs(float(c)) for c in (P + Q)))
    info['gap_M'] = float(abs(dl - arc.r) / M)
    if h2 < 0:
        return [], info
    L = mp.sqrt(mpf(L2))
    hl = mp.sqrt(mpf(h2) / mpf(L2))      # h / L in parameter units
    cands = [(mpf(s) + hl), (mpf(s) - hl)] if h2 != 0 else [mpf(s)]
    pts = []
    margins = []
    for t in cands:
        X = (mpf(P[0]) + t * mpf(dx), mpf(P[1]) + t * mpf(dy))
        snapped = None
        for E, te in ((P, 0), (Q, 1)):
            if arc.on_circle_exact(E) and abs(t - te) < mp.mpf(10) ** -40:
                snapped = E
        for E in (arc.A, arc.C, arc.B):
            # arc control point exactly on the segment line
            if cross(P, Q, E) == 0 and arc.on_circle_exact(E):
                Em = (mpf(E[0]), mpf(E[1]))
                if mp.sqrt((Em[0] - X[0]) ** 2 + (Em[1] - X[1]) ** 2) <= arc.r * mp.mpf(10) ** -40:
                    snapped = E
        if snapped is not None:
            tt = ((snapped[0] - P[0]) * dx + (snapped[1] - P[1]) * dy) / L2
            if 0 <= tt <= 1 and arc.contains_exact(snapped):
                pts.append((mpf(snapped[0]), mpf(snapped[1])))
            continue
        on_seg = (t >= 0 and t <= 1)
        if on_seg and arc.in_span_mp(X):
            pts.append(X)
        margins.append(min(arc.span_margin(X), min(abs(t), abs(1 - t)) * L / arc.r))
    info['margin'] = float(min(margins)) if margins else float('inf')
    return pts, info


# ---------------------------------------------------------------------------
# Point-to-arc distance (theories/LECArcRow.v arc_dist, arc_dist_exact)
# ---------------------------------------------------------------------------

def arc_point_distance(arc, P):
    Pm = (mpf(Fr(P[0])), mpf(Fr(P[1])))
    dA = mp.sqrt((Pm[0] - mpf(arc.A[0])) ** 2 + (Pm[1] - mpf(arc.A[1])) ** 2)
    dC = mp.sqrt((Pm[0] - mpf(arc.C[0])) ** 2 + (Pm[1] - mpf(arc.C[1])) ** 2)
    ends = min(dA, dC)
    dO = mp.sqrt((Pm[0] - arc.Om[0]) ** 2 + (Pm[1] - arc.Om[1]) ** 2)
    if dO == 0:
        return arc.r                          # centre case: every arc point at r
    foot = (arc.Om[0] + arc.r * (Pm[0] - arc.Om[0]) / dO,
            arc.Om[1] + arc.r * (Pm[1] - arc.Om[1]) / dO)
    if arc.in_span_mp(foot):
        return min(abs(dO - arc.r), ends)
    return ends


def point_on_arc_dist(arc, X):
    """How far a computed point X (floats) is from the arc (mp)."""
    return arc_point_distance(arc, X)


def seg_point_distance(P, Q, X):
    P = (mpf(Fr(P[0])), mpf(Fr(P[1]))); Q = (mpf(Fr(Q[0])), mpf(Fr(Q[1])))
    X = (mpf(Fr(X[0])), mpf(Fr(X[1])))
    dx, dy = Q[0] - P[0], Q[1] - P[1]
    L2 = dx * dx + dy * dy
    t = 0 if L2 == 0 else max(0, min(1, ((X[0] - P[0]) * dx + (X[1] - P[1]) * dy) / L2))
    fx, fy = P[0] + t * dx, P[1] + t * dy
    return mp.sqrt((X[0] - fx) ** 2 + (X[1] - fy) ** 2)
