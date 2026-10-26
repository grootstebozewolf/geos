/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2026 Jeroen Bloemscheer
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <algorithm/ClothoidGeometry.h>
#include <algorithm/ClothoidMath.h>

#include <geos/constants.h>
#include <geos/geom/CircularArc.h>
#include <geos/geom/Clothoid.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/LineSegment.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

using geos::geom::CircularArc;
using geos::geom::Clothoid;
using geos::geom::CoordinateXY;
using geos::geom::Envelope;
using geos::geom::LineSegment;

namespace geos {
namespace algorithm {
namespace clothoid {

namespace {

// Endpoint snap tolerance of Clothoid, in ulps of the reference location extent.
constexpr double TOLERANCE_ULPS = 32.0;

// Past this many halvings a span is visited as it is.
constexpr int MAX_DEPTH = 64;

// Clothoid against clothoid refines each candidate span once it is this
// fraction of the curve's extent.
constexpr double PAIR_LEAF_FRACTION = 1.0 / 1024.0;

double
dot(const CoordinateXY& a, const CoordinateXY& b)
{
    return a.x * b.x + a.y * b.y;
}

CoordinateXY
minus(const CoordinateXY& a, const CoordinateXY& b)
{
    return CoordinateXY(a.x - b.x, a.y - b.y);
}

/// Distance from q to the box spanned by a and b.
double
boxDistance(const CoordinateXY& q, const CoordinateXY& a, const CoordinateXY& b)
{
    const double dx = std::max(0.0, std::max(std::min(a.x, b.x) - q.x, q.x - std::max(a.x, b.x)));
    const double dy = std::max(0.0, std::max(std::min(a.y, b.y) - q.y, q.y - std::max(a.y, b.y)));
    return std::hypot(dx, dy);
}

/// Distance from q to the farthest corner of the box spanned by a and b.
double
boxFarthest(const CoordinateXY& q, const CoordinateXY& a, const CoordinateXY& b)
{
    const double dx = std::max(std::abs(q.x - a.x), std::abs(q.x - b.x));
    const double dy = std::max(std::abs(q.y - a.y), std::abs(q.y - b.y));
    return std::hypot(dx, dy);
}

/// Bisects a sign change of f on [a, b], where f(a) = fa.
template<typename F>
double
bisect(F&& f, double a, double b, double fa)
{
    for (int i = 0; i < 2 * MAX_DEPTH; i++) {
        const double m = a + 0.5 * (b - a);
        if (m == a || m == b) {
            break;
        }
        const double fm = f(m);
        if (fm == 0.0) {
            return m;
        }
        if ((fm > 0.0) == (fa > 0.0)) {
            a = m;
            fa = fm;
        }
        else {
            b = m;
        }
    }
    return a + 0.5 * (b - a);
}

bool
signChange(double fa, double fb)
{
    return (fa > 0.0 && fb < 0.0) || (fa < 0.0 && fb > 0.0);
}

struct Span {
    double a;
    double b;
    CoordinateXY pa;
    CoordinateXY pb;
};

/// Splits the curve into spans on which |P(t) - c| is monotone, and visits
/// them. The critical points of |P - c| are the roots of g = (P - c).T; a
/// span is accepted when a bound on g'' shows g is monotone on it, and split
/// at its root when g changes sign. Spans for which prune(pa, pb) holds are
/// skipped. A span is a sub-span of a monotone piece, so the box of its
/// endpoints holds it.
template<typename Prune, typename Visit>
class CriticalWalk {
public:
    CriticalWalk(const PlacedClothoid& curve, const CoordinateXY& c, Prune& prune, Visit& visit)
        : m_curve(curve), m_c(c), m_prune(prune), m_visit(visit)
    {
        m_invA2 = std::abs(curve.headingRate(1.0));
        m_sigma = curve.maxSingularValue();
    }

    void run()
    {
        const auto& breaks = m_curve.monotoneBreaks();
        for (std::size_t i = 1; i < breaks.size(); i++) {
            const double a = breaks[i - 1];
            const double b = breaks[i];
            const CoordinateXY pa = m_curve.pointAt(a);
            const CoordinateXY pb = m_curve.pointAt(b);
            walk(Span{a, b, pa, pb}, g(a, pa), g(b, pb), 0);
        }
    }

private:
    double g(double t, const CoordinateXY& p) const
    {
        return dot(minus(p, m_c), m_curve.tangentAt(t));
    }

    void walk(const Span& s, double ga, double gb, int depth)
    {
        if (m_prune(s.pa, s.pb)) {
            return;
        }
        const double m = s.a + 0.5 * (s.b - s.a);
        if (m == s.a || m == s.b) {
            m_visit(s);
            return;
        }
        const CoordinateXY pm = m_curve.pointAt(m);
        const CoordinateXY tm = m_curve.tangentAt(m);
        const CoordinateXY nm = m_curve.normalAt(m);
        const CoordinateXY rm = minus(pm, m_c);
        const double gm = dot(rm, tm);
        // g' = |T|^2 + theta' (P - c).N and
        // g'' = 3 theta' T.N + theta'' (P - c).N - theta'^2 (P - c).T.
        const double slope = dot(tm, tm) + m_curve.headingRate(m) * dot(rm, nm);
        const double rate = std::max(std::abs(s.a), std::abs(s.b)) * m_invA2;
        const double reach = boxFarthest(m_c, s.pa, s.pb);
        const double bound = 3.0 * rate * m_sigma * m_sigma
                             + m_invA2 * reach * m_sigma
                             + rate * rate * reach * m_sigma;
        if (std::abs(slope) > 0.5 * std::abs(s.b - s.a) * bound) {
            // g is monotone on the span: at most one critical point.
            if (signChange(ga, gb)) {
                const double r = bisect([&](double t) {
                    return g(t, m_curve.pointAt(t));
                }, s.a, s.b, ga);
                const CoordinateXY pr = m_curve.pointAt(r);
                m_visit(Span{s.a, r, s.pa, pr});
                m_visit(Span{r, s.b, pr, s.pb});
            }
            else {
                m_visit(s);
            }
            return;
        }
        if (depth >= MAX_DEPTH) {
            m_visit(Span{s.a, m, s.pa, pm});
            m_visit(Span{m, s.b, pm, s.pb});
            return;
        }
        const Span left{s.a, m, s.pa, pm};
        const Span right{m, s.b, pm, s.pb};
        if (boxDistance(m_c, s.pa, pm) <= boxDistance(m_c, pm, s.pb)) {
            walk(left, ga, gm, depth + 1);
            walk(right, gm, gb, depth + 1);
        }
        else {
            walk(right, gm, gb, depth + 1);
            walk(left, ga, gm, depth + 1);
        }
    }

    const PlacedClothoid& m_curve;
    CoordinateXY m_c;
    Prune& m_prune;
    Visit& m_visit;
    double m_invA2;
    double m_sigma;
};

template<typename Prune, typename Visit>
void
walkCritical(const PlacedClothoid& curve, const CoordinateXY& c, Prune& prune, Visit& visit)
{
    CriticalWalk<Prune, Visit> walk(curve, c, prune, visit);
    walk.run();
}

class Best {
public:
    Best()
    {
        m_nearest.distance = std::numeric_limits<double>::infinity();
    }

    void consider(double distance, const CoordinateXY& onCurve, const CoordinateXY& onTarget)
    {
        if (distance < m_nearest.distance) {
            m_nearest.distance = distance;
            m_nearest.onCurve = onCurve;
            m_nearest.onTarget = onTarget;
        }
    }

    void consider(const Nearest& n)
    {
        consider(n.distance, n.onCurve, n.onTarget);
    }

    double distance() const
    {
        return m_nearest.distance;
    }

    const Nearest& nearest() const
    {
        return m_nearest;
    }

private:
    Nearest m_nearest;
};

Nearest
swapped(const Nearest& n)
{
    Nearest s;
    s.distance = n.distance;
    s.onCurve = n.onTarget;
    s.onTarget = n.onCurve;
    return s;
}

/// Report a pair within the evaluation tolerance as coincident at the curve
/// point, as PointLocator treats such points as on the curve.
Nearest
snapped(const Nearest& n, double tolerance)
{
    if (n.distance > tolerance || n.distance == 0.0) {
        return n;
    }
    Nearest s;
    s.distance = 0.0;
    s.onCurve = n.onCurve;
    s.onTarget = n.onCurve;
    return s;
}

/// Distance between the boxes spanned by (a, b) and env.
double
boxEnvelopeDistance(const CoordinateXY& a, const CoordinateXY& b, const Envelope& env)
{
    return Envelope(a, b).distance(env);
}

} // namespace

PlacedClothoid::PlacedClothoid(const Clothoid& curve)
{
    const geom::ClothoidReferenceLocation loc = curve.getReferenceLocation();
    m_ox = loc.origin.x;
    m_oy = loc.origin.y;
    m_ux = loc.u.x;
    m_uy = loc.u.y;
    m_vx = loc.v.x;
    m_vy = loc.v.y;
    m_scale = std::abs(curve.getScaleFactor());
    m_t0 = curve.getStartDistance();
    m_t1 = curve.getEndDistance();
    m_start = curve.getCoordinatesRO()->getAt<CoordinateXY>(0);
    m_end = curve.getCoordinatesRO()->getAt<CoordinateXY>(1);
    m_envelope = *curve.getEnvelopeInternal();

    const double a = m_ux * m_ux + m_uy * m_uy;
    const double b = m_vx * m_vx + m_vy * m_vy;
    const double c = m_ux * m_vx + m_uy * m_vy;
    const double disc = std::sqrt(std::max(0.0, (a - b) * (a - b) + 4.0 * c * c));
    m_sigmaMax = std::sqrt(0.5 * (a + b + disc));

    const double extent = std::max(std::abs(m_ox), std::abs(m_oy))
                          + m_sigmaMax * std::max(std::abs(m_t0), std::abs(m_t1));
    m_tolerance = TOLERANCE_ULPS * std::numeric_limits<double>::epsilon() * extent;

    m_breaks = perpendicularTangents(1.0, 0.0);
    const std::vector<double> yBreaks = perpendicularTangents(0.0, 1.0);
    m_breaks.insert(m_breaks.end(), yBreaks.begin(), yBreaks.end());
    m_breaks.push_back(std::min(m_t0, m_t1));
    m_breaks.push_back(std::max(m_t0, m_t1));
    std::sort(m_breaks.begin(), m_breaks.end());
    m_breaks.erase(std::unique(m_breaks.begin(), m_breaks.end()), m_breaks.end());
    if (m_t1 < m_t0) {
        std::reverse(m_breaks.begin(), m_breaks.end());
    }
}

CoordinateXY
PlacedClothoid::pointAt(double t) const
{
    if (t == m_t0) {
        return m_start;
    }
    if (t == m_t1) {
        return m_end;
    }
    const CoordinateXY c = constructiveOffset(t, m_scale);
    return CoordinateXY(m_ox + c.x * m_ux + c.y * m_vx,
                        m_oy + c.x * m_uy + c.y * m_vy);
}

CoordinateXY
PlacedClothoid::tangentAt(double t) const
{
    double c = 0.0;
    double s = 0.0;
    halfSquarePhase(t, m_scale, c, s);
    return CoordinateXY(m_ux * c + m_vx * s, m_uy * c + m_vy * s);
}

CoordinateXY
PlacedClothoid::normalAt(double t) const
{
    double c = 0.0;
    double s = 0.0;
    halfSquarePhase(t, m_scale, c, s);
    return CoordinateXY(-m_ux * s + m_vx * c, -m_uy * s + m_vy * c);
}

double
PlacedClothoid::headingRate(double t) const
{
    return t / (m_scale * m_scale);
}

std::vector<double>
PlacedClothoid::perpendicularTangents(double dx, double dy) const
{
    // T.d = (u.d) cos theta + (v.d) sin theta vanishes where theta is
    // base + k pi. Each side of the inflection is walked from its own start
    // with the reduced heading there, as the envelope does.
    std::vector<double> out;
    const double coeffCos = m_ux * dx + m_uy * dy;
    const double coeffSin = m_vx * dx + m_vy * dy;
    if (!(std::hypot(coeffCos, coeffSin) > 0.0)) {
        return out;
    }
    const double base = std::atan2(coeffSin, coeffCos) + MATH_PI / 2.0;
    const double tLo = std::min(m_t0, m_t1);
    const double tHi = std::max(m_t0, m_t1);
    const double a2 = m_scale * m_scale;

    auto solveSide = [&](double s0, double s1, double sign) {
        if (!(s1 > s0)) {
            return;
        }
        double c = 0.0;
        double s = 0.0;
        halfSquarePhase(s0, m_scale, c, s);
        const double turn = (s1 - s0) * (s1 + s0) / (2.0 * a2);
        double offset = std::fmod(base - std::atan2(s, c), MATH_PI);
        if (offset < 0.0) {
            offset += MATH_PI;
        }
        for (double k = 0.0; ; k += 1.0) {
            const double delta = offset + k * MATH_PI;
            if (!(delta <= turn)) {
                break;
            }
            const double w = 2.0 * a2 * delta;
            const double t = sign * ((w > 0.0) ? s0 + w / (s0 + std::sqrt(s0 * s0 + w)) : s0);
            if (t > tLo && t < tHi) {
                out.push_back(t);
            }
        }
    };
    if (tHi > 0.0) {
        solveSide(std::max(0.0, tLo), tHi, 1.0);
    }
    if (tLo < 0.0) {
        solveSide(std::max(0.0, -tHi), -tLo, -1.0);
    }
    std::sort(out.begin(), out.end());
    return out;
}

Nearest
nearestToPoint(const PlacedClothoid& curve, const CoordinateXY& q)
{
    Best best;
    const CoordinateXY start = curve.pointAt(curve.startDistance());
    const CoordinateXY end = curve.pointAt(curve.endDistance());
    best.consider(start.distance(q), start, q);
    best.consider(end.distance(q), end, q);

    // |P - q| is monotone on each visited span, so its minimum there is at an end.
    auto prune = [&](const CoordinateXY& pa, const CoordinateXY& pb) {
        return boxDistance(q, pa, pb) >= best.distance();
    };
    auto visit = [&](const Span& s) {
        best.consider(s.pa.distance(q), s.pa, q);
        best.consider(s.pb.distance(q), s.pb, q);
    };
    walkCritical(curve, q, prune, visit);
    return snapped(best.nearest(), curve.tolerance());
}

Nearest
nearestToSegment(const PlacedClothoid& curve, const CoordinateXY& p0, const CoordinateXY& p1)
{
    if (p0.equals2D(p1)) {
        return nearestToPoint(curve, p0);
    }
    Best best;
    best.consider(nearestToPoint(curve, p0));
    best.consider(nearestToPoint(curve, p1));

    const double tStart = curve.startDistance();
    const double tEnd = curve.endDistance();
    for (double t : {tStart, tEnd}) {
        const CoordinateXY p = curve.pointAt(t);
        const CoordinateXY onSegment = LineSegment::closestPoint(p0, p1, p);
        best.consider(p.distance(onSegment), p, onSegment);
    }

    // nu = (P - p0).n is monotone between the points where the tangent is
    // parallel to the segment. Its extremes there are the interior
    // candidates, and a sign change inside the segment is a crossing.
    const CoordinateXY d = minus(p1, p0);
    const double length2 = dot(d, d);
    const double length = std::sqrt(length2);
    const CoordinateXY n(-d.y / length, d.x / length);
    auto nu = [&](const CoordinateXY& p) {
        return dot(minus(p, p0), n);
    };
    auto along = [&](const CoordinateXY& p) {
        return dot(minus(p, p0), d) / length2;
    };

    std::vector<double> stops = curve.perpendicularTangents(n.x, n.y);
    for (double t : stops) {
        const CoordinateXY p = curve.pointAt(t);
        const double s = along(p);
        if (s >= 0.0 && s <= 1.0) {
            const CoordinateXY onSegment(p0.x + s * d.x, p0.y + s * d.y);
            best.consider(std::abs(nu(p)), p, onSegment);
        }
    }
    stops.insert(stops.begin(), std::min(tStart, tEnd));
    stops.push_back(std::max(tStart, tEnd));
    for (std::size_t i = 1; i < stops.size() && best.distance() > 0.0; i++) {
        const double a = stops[i - 1];
        const double b = stops[i];
        const double nuA = nu(curve.pointAt(a));
        if (!signChange(nuA, nu(curve.pointAt(b)))) {
            continue;
        }
        const double t = bisect([&](double x) {
            return nu(curve.pointAt(x));
        }, a, b, nuA);
        const CoordinateXY p = curve.pointAt(t);
        const double s = along(p);
        if (s >= 0.0 && s <= 1.0) {
            best.consider(0.0, p, p);
        }
    }
    return snapped(best.nearest(), curve.tolerance());
}

Nearest
nearestToArc(const PlacedClothoid& curve, const CircularArc& arc)
{
    if (arc.isLinear()) {
        Best best;
        best.consider(nearestToSegment(curve, arc.p0(), arc.p1()));
        best.consider(nearestToSegment(curve, arc.p1(), arc.p2()));
        return best.nearest();
    }

    Best best;
    best.consider(nearestToPoint(curve, arc.p0()));
    best.consider(nearestToPoint(curve, arc.p2()));

    const CoordinateXY& c = arc.getCenter();
    const double r = arc.getRadius();
    const Envelope arcEnv = arc.getEnvelope();

    auto considerPoint = [&](const CoordinateXY& p) {
        const CoordinateXY onArc = arc.closestPoint(p);
        best.consider(p.distance(onArc), p, onArc);
    };
    // Within a span |P - c| is monotone, so h = |P - c| - r crosses zero at
    // most once; elsewhere the distance to the arc is that to its circle,
    // at a span end, or to an arc endpoint.
    auto prune = [&](const CoordinateXY& pa, const CoordinateXY& pb) {
        const double nearC = boxDistance(c, pa, pb);
        const double farC = boxFarthest(c, pa, pb);
        const double circle = std::max(0.0, std::max(nearC - r, r - farC));
        const double box = boxEnvelopeDistance(pa, pb, arcEnv);
        return std::max(circle, box) >= best.distance();
    };
    auto h = [&](const CoordinateXY& p) {
        return p.distance(c) - r;
    };
    auto visit = [&](const Span& s) {
        considerPoint(s.pa);
        considerPoint(s.pb);
        const double ha = h(s.pa);
        if (best.distance() > 0.0 && signChange(ha, h(s.pb))) {
            const double t = bisect([&](double x) {
                return h(curve.pointAt(x));
            }, s.a, s.b, ha);
            const CoordinateXY p = curve.pointAt(t);
            if (arc.containsPointOnCircle(p)) {
                best.consider(0.0, p, p);
            }
        }
    };
    walkCritical(curve, c, prune, visit);
    return snapped(best.nearest(), curve.tolerance());
}

Nearest
nearestToClothoid(const PlacedClothoid& curve, const PlacedClothoid& other)
{
    Best best;
    best.consider(swapped(nearestToPoint(other, curve.pointAt(curve.startDistance()))));
    best.consider(swapped(nearestToPoint(other, curve.pointAt(curve.endDistance()))));
    best.consider(nearestToPoint(curve, other.pointAt(other.startDistance())));
    best.consider(nearestToPoint(curve, other.pointAt(other.endDistance())));

    // F(t) = distance from P(t) to the other curve is 1-Lipschitz in arc
    // length, and a monotone span is no longer than the L1 distance between
    // its ends. Spans that cannot beat the best are dropped, and the rest
    // are refined by golden section once they are small.
    struct Item {
        double lower;
        double a;
        double b;
        CoordinateXY pa;
        CoordinateXY pb;
        double fa;
        double fb;
        bool operator<(const Item& o) const
        {
            return lower > o.lower;
        }
    };
    auto evaluate = [&](double t) {
        const CoordinateXY p = curve.pointAt(t);
        const Nearest n = swapped(nearestToPoint(other, p));
        best.consider(n);
        return n.distance;
    };
    auto lowerBound = [&](const CoordinateXY& pa, const CoordinateXY& pb, double fa, double fb) {
        const double reach = std::abs(pb.x - pa.x) + std::abs(pb.y - pa.y);
        return std::max(0.5 * (fa + fb - reach), boxEnvelopeDistance(pa, pb, other.envelope()));
    };

    const Envelope& env = curve.envelope();
    const double leaf = PAIR_LEAF_FRACTION * std::hypot(env.getWidth(), env.getHeight());

    std::priority_queue<Item> queue;
    const auto& breaks = curve.monotoneBreaks();
    for (std::size_t i = 1; i < breaks.size(); i++) {
        const double a = breaks[i - 1];
        const double b = breaks[i];
        const CoordinateXY pa = curve.pointAt(a);
        const CoordinateXY pb = curve.pointAt(b);
        if (boxEnvelopeDistance(pa, pb, other.envelope()) >= best.distance()) {
            continue;
        }
        const double fa = evaluate(a);
        const double fb = evaluate(b);
        queue.push(Item{lowerBound(pa, pb, fa, fb), a, b, pa, pb, fa, fb});
    }

    const double invPhi = 0.5 * (std::sqrt(5.0) - 1.0);
    while (!queue.empty() && best.distance() > 0.0) {
        const Item item = queue.top();
        queue.pop();
        if (item.lower >= best.distance()) {
            break;
        }
        const double reach = std::abs(item.pb.x - item.pa.x) + std::abs(item.pb.y - item.pa.y);
        if (reach <= leaf) {
            double lo = item.a;
            double hi = item.b;
            double x1 = hi - invPhi * (hi - lo);
            double x2 = lo + invPhi * (hi - lo);
            double f1 = evaluate(x1);
            double f2 = evaluate(x2);
            for (int i = 0; i < 2 * MAX_DEPTH && x1 != x2; i++) {
                if (f1 <= f2) {
                    hi = x2;
                    x2 = x1;
                    f2 = f1;
                    x1 = hi - invPhi * (hi - lo);
                    f1 = evaluate(x1);
                }
                else {
                    lo = x1;
                    x1 = x2;
                    f1 = f2;
                    x2 = lo + invPhi * (hi - lo);
                    f2 = evaluate(x2);
                }
            }
            continue;
        }
        const double m = item.a + 0.5 * (item.b - item.a);
        if (m == item.a || m == item.b) {
            continue;
        }
        const CoordinateXY pm = curve.pointAt(m);
        const double fm = evaluate(m);
        queue.push(Item{lowerBound(item.pa, pm, item.fa, fm), item.a, m, item.pa, pm, item.fa, fm});
        queue.push(Item{lowerBound(pm, item.pb, fm, item.fb), m, item.b, pm, item.pb, fm, item.fb});
    }
    return snapped(best.nearest(), std::max(curve.tolerance(), other.tolerance()));
}

void
countRayCrossings(const PlacedClothoid& curve, const CoordinateXY& q,
                  std::size_t& crossings, bool& onCurve)
{
    const Envelope& env = curve.envelope();
    const double tol = curve.tolerance();
    if (env.getMaxX() < q.x - tol || env.getMaxY() < q.y - tol || env.getMinY() > q.y + tol) {
        return;
    }
    if (env.getMinX() <= q.x + tol && nearestToPoint(curve, q).distance <= tol) {
        onCurve = true;
        return;
    }

    // Each monotone piece crosses the horizontal line through q at most once.
    // An upward piece includes its start and excludes its end; a downward
    // piece excludes its start and includes its end.
    const auto& breaks = curve.monotoneBreaks();
    for (std::size_t i = 1; i < breaks.size(); i++) {
        const double a = breaks[i - 1];
        const double b = breaks[i];
        const CoordinateXY pa = curve.pointAt(a);
        const CoordinateXY pb = curve.pointAt(b);
        if (pa.x < q.x && pb.x < q.x) {
            continue;
        }
        if (!((pa.y > q.y && pb.y <= q.y) || (pb.y > q.y && pa.y <= q.y))) {
            continue;
        }
        double x;
        if (pa.y == q.y) {
            x = pa.x;
        }
        else if (pb.y == q.y) {
            x = pb.x;
        }
        else {
            const double t = bisect([&](double s) {
                return curve.pointAt(s).y - q.y;
            }, a, b, pa.y - q.y);
            x = curve.pointAt(t).x;
        }
        if (x > q.x) {
            crossings++;
        }
    }
}

}
}
}
