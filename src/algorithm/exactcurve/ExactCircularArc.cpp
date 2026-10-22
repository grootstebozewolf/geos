/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/algorithm/exactcurve/ExactCircularArc.h>
#include <geos/algorithm/CircularArcs.h>
#include <geos/algorithm/Orientation.h>
#include <geos/constants.h>
#include <geos/geom/CircularString.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/math/DD.h>
#include <geos/util/IllegalArgumentException.h>

#include <algorithm>
#include <cmath>
#include <limits>

using geos::geom::CircularString;
using geos::geom::CoordinateSequence;
using geos::geom::CoordinateXY;
using geos::geom::CoordinateXYZM;
using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::math::DD;

namespace geos {
namespace algorithm {
namespace exactcurve {

namespace {

constexpr double TWO_PI = 2.0 * MATH_PI;
constexpr double DEFAULT_TOLERANCE_FRACTION = 0.01;

double
normPos(double angle)
{
    angle = std::fmod(angle, TWO_PI);
    if (angle < 0.0) {
        angle += TWO_PI;
    }
    return angle;
}

double
signedShort(double ux, double uy, double vx, double vy)
{
    return std::atan2(ux * vy - uy * vx, ux * vx + uy * vy);
}

} // namespace

ExactCircularArc::ExactCircularArc(const CoordinateXY& start,
                                   const CoordinateXY& mid,
                                   const CoordinateXY& end)
    : m_start(start)
    , m_mid(mid)
    , m_end(end)
    , m_cx(std::numeric_limits<double>::quiet_NaN())
    , m_cy(std::numeric_limits<double>::quiet_NaN())
    , m_r(0.0)
    , m_a0(0.0)
    , m_ccw(true)
    , m_sweep(0.0)
    , m_arc(false)
{
    if (!tryCircumcircle(start, mid, end, m_cx, m_cy, m_r)) {
        return;
    }
    m_a0 = std::atan2(start.y - m_cy, start.x - m_cx);
    const double shortSE = signedShort(start.x - m_cx, start.y - m_cy,
                                       end.x - m_cx, end.y - m_cy);
    const double shortSM = signedShort(start.x - m_cx, start.y - m_cy,
                                       mid.x - m_cx, mid.y - m_cy);
    m_ccw = normPos(shortSM) < normPos(shortSE);
    m_sweep = m_ccw ? normPos(shortSE) : normPos(-shortSE);
    if (m_sweep == 0.0) {
        m_sweep = TWO_PI;
    }
    m_arc = true;
}

double
ExactCircularArc::lengthOf(const CoordinateXY& start,
                           const CoordinateXY& mid,
                           const CoordinateXY& end)
{
    double cx, cy, r;
    if (!tryCircumcircle(start, mid, end, cx, cy, r)) {
        return start.distance(end);
    }
    return r * directedSweep(cx, cy, start, mid, end);
}

// Port of JTS 9797c2c4.
double
ExactCircularArc::length() const
{
    return m_arc ? m_r * m_sweep : m_start.distance(m_end);
}

CoordinateXY
ExactCircularArc::center() const
{
    return m_arc ? CoordinateXY{m_cx, m_cy} : CoordinateXY{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
}

CoordinateXY
ExactCircularArc::pointAt(double t) const
{
    if (!std::isfinite(t) || t < 0.0 || t > 1.0) {
        throw util::IllegalArgumentException("t must be in [0,1]");
    }
    if (t == 0.0) {
        return m_start;
    }
    if (t == 1.0) {
        return m_end;
    }
    if (!m_arc) {
        return CoordinateXY{
            m_start.x + t * (m_end.x - m_start.x),
            m_start.y + t * (m_end.y - m_start.y)};
    }
    const double ang = m_a0 + (m_ccw ? m_sweep : -m_sweep) * t;
    return CoordinateXY{m_cx + m_r * std::cos(ang), m_cy + m_r * std::sin(ang)};
}

std::unique_ptr<Geometry>
ExactCircularArc::toLinear(double tolerance) const
{
    if (tolerance < 0.0) {
        throw util::IllegalArgumentException("tolerance must be non-negative");
    }
    const GeometryFactory* gf = GeometryFactory::getDefaultInstance();
    if (!m_arc) {
        auto seq = std::make_unique<CoordinateSequence>(2, false, false);
        seq->setAt(m_start, 0);
        seq->setAt(m_end, 1);
        return gf->createLineString(std::move(seq));
    }
    const double eps = tolerance == 0.0 ? m_r * DEFAULT_TOLERANCE_FRACTION : tolerance;
    const int segments = segmentCount(m_r, m_sweep, eps);
    double delta = m_sweep / segments;
    if (!m_ccw) {
        delta = -delta;
    }
    auto seq = std::make_unique<CoordinateSequence>(
        static_cast<std::size_t>(segments + 1), false, false);
    seq->setAt(m_start, 0);
    for (int i = 1; i < segments; ++i) {
        const double ang = m_a0 + i * delta;
        seq->setAt(CoordinateXY{m_cx + m_r * std::cos(ang),
                                m_cy + m_r * std::sin(ang)},
                   static_cast<std::size_t>(i));
    }
    seq->setAt(m_end, static_cast<std::size_t>(segments));
    return gf->createLineString(std::move(seq));
}

bool
ExactCircularArc::chordLeArc() const
{
    const double chord = chordLength();
    if (!m_arc) {
        return true;
    }
    const double arcLen = m_r * m_sweep;
    if (chord <= arcLen) {
        return true;
    }
    const double chordFromSweep = 2.0 * m_r * std::sin(0.5 * m_sweep);
    const double bound = std::max(arcLen, chordFromSweep);
    return chord <= bound + ulp(std::max(bound, chord));
}

bool
ExactCircularArc::inArc(const CoordinateXY& p, double radialTol) const
{
    if (!m_arc) {
        return onSegment(p, m_start, m_end, radialTol);
    }
    const double dx = p.x - m_cx;
    const double dy = p.y - m_cy;
    const double d2 = dx * dx + dy * dy;
    const double r2 = m_r * m_r;
    const double tol2 = radialTol * (2.0 * m_r + radialTol);
    if (std::fabs(d2 - r2) > tol2) {
        return false;
    }
    return onSweep(p);
}

double
ExactCircularArc::circularSegmentArea() const
{
    if (!m_arc) {
        return 0.0;
    }
    return 0.5 * m_r * m_r * (m_sweep - std::sin(m_sweep));
}

CoordinateXY
ExactCircularArc::arcLengthCentroid() const
{
    if (!m_arc) {
        return CoordinateXY{0.5 * (m_start.x + m_end.x),
                            0.5 * (m_start.y + m_end.y)};
    }
    if (m_sweep == 0.0) {
        return m_start;
    }
    const double signedSweep = m_ccw ? m_sweep : -m_sweep;
    const double a1 = m_a0 + signedSweep;
    const double k = m_r / signedSweep;
    return CoordinateXY{
        m_cx + k * (std::sin(a1) - std::sin(m_a0)),
        m_cy + k * (-std::cos(a1) + std::cos(m_a0))};
}

bool
ExactCircularArc::tryCircumcircle(const CoordinateXY& a,
                                  const CoordinateXY& b,
                                  const CoordinateXY& c,
                                  double& cx, double& cy, double& r)
{
    if (Orientation::index(a, b, c) == Orientation::COLLINEAR) {
        return false;
    }
    const double d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if (d == 0.0) {
        return false;
    }
    const double a2 = a.x * a.x + a.y * a.y;
    const double b2 = b.x * b.x + b.y * b.y;
    const double c2 = c.x * c.x + c.y * c.y;
    cx = (a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / d;
    cy = (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / d;
    r = std::hypot(a.x - cx, a.y - cy);
    return std::isfinite(r) && r != 0.0;
}

int
ExactCircularArc::segmentCount(double radius, double sweep, double tolerance)
{
    if (tolerance >= radius) {
        return 1;
    }
    const double thetaMax = 2.0 * std::acos(1.0 - tolerance / radius);
    if (!std::isfinite(thetaMax) || thetaMax <= 0.0) {
        return 1;
    }
    const int n = static_cast<int>(std::ceil(sweep / thetaMax));
    return n < 1 ? 1 : n;
}

double
ExactCircularArc::directedSweep(double cx, double cy,
                                const CoordinateXY& start,
                                const CoordinateXY& mid,
                                const CoordinateXY& end)
{
    const double shortSE = signedShort(start.x - cx, start.y - cy,
                                       end.x - cx, end.y - cy);
    const double shortSM = signedShort(start.x - cx, start.y - cy,
                                       mid.x - cx, mid.y - cy);
    const bool ccw = normPos(shortSM) < normPos(shortSE);
    double sweep = ccw ? normPos(shortSE) : normPos(-shortSE);
    return sweep == 0.0 ? TWO_PI : sweep;
}

bool
ExactCircularArc::onSweep(const CoordinateXY& p) const
{
    if (!m_arc) {
        return false;
    }
    const double s = signedShort(m_start.x - m_cx, m_start.y - m_cy,
                                 p.x - m_cx, p.y - m_cy);
    const double travelled = m_ccw ? normPos(s) : normPos(-s);
    return travelled <= m_sweep + ulp(m_sweep);
}

bool
ExactCircularArc::onSegment(const CoordinateXY& p,
                            const CoordinateXY& a,
                            const CoordinateXY& b,
                            double tol)
{
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double len2 = dx * dx + dy * dy;
    if (len2 == 0.0) {
        return p.distance(a) <= tol;
    }
    double t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
    if (t < 0.0) {
        t = 0.0;
    }
    else if (t > 1.0) {
        t = 1.0;
    }
    const double px = a.x + t * dx - p.x;
    const double py = a.y + t * dy - p.y;
    return px * px + py * py <= tol * tol;
}

double
ExactCircularArc::ulp(double value)
{
    if (!std::isfinite(value)) {
        return value;
    }
    const double abs = std::fabs(value);
    return std::nextafter(abs, std::numeric_limits<double>::infinity()) - abs;
}

bool
ExactCircularArc::isProperInterior(const CoordinateXY& p, double radialTol) const
{
    if (!m_arc) {
        return false;
    }
    if (!inArc(p, radialTol)) {
        return false;
    }
    const double endTol = std::max(1.0e-12, 1.0e-12 * m_r);
    if (p.distance(m_start) <= endTol || p.distance(m_end) <= endTol) {
        return false;
    }
    return true;
}

CoordinateXY
ExactCircularArc::angularMid(const CoordinateXY& from,
                             const CoordinateXY& to) const
{
    return CircularArcs::getMidpoint(from, to, center(), m_r, m_ccw);
}

std::optional<Hit>
ExactCircularArc::hit(const ExactCircularArc& a, const ExactCircularArc& b)
{
    if (!a.isArc() || !b.isArc()) {
        return std::nullopt;
    }

    const bool swapArgs = a.center().compareTo(b.center()) > 0;
    const ExactCircularArc& lhs = swapArgs ? b : a;
    const ExactCircularArc& rhs = swapArgs ? a : b;

    auto [c1x, c1y] = CircularArcs::getCenterDD(lhs.getStart(), lhs.getMid(), lhs.getEnd());
    auto [c2x, c2y] = CircularArcs::getCenterDD(rhs.getStart(), rhs.getMid(), rhs.getEnd());
    const DD r1x = c1x - lhs.getStart().x;
    const DD r1y = c1y - lhs.getStart().y;
    const DD r2x = c2x - rhs.getStart().x;
    const DD r2y = c2y - rhs.getStart().y;
    const DD r1sq = r1x * r1x + r1y * r1y;
    const DD r2sq = r2x * r2x + r2y * r2y;

    CoordinateXY isect0, isect1;
    const int n = CircularArcs::circleIntersectsCircle(
        c1x, c1y, r1sq, c2x, c2y, r2sq, isect0, isect1);
    if (n <= 0) {
        // n < 0 is cocircular: not a proper crossing.
        return std::nullopt;
    }

    const double tol = 1.0e-8 * std::max(1.0, std::max(a.radius(), b.radius()));
    CoordinateXY found;
    int nProper = 0;
    const CoordinateXY cands[2] = {isect0, isect1};
    for (int i = 0; i < n; ++i) {
        if (a.isProperInterior(cands[i], tol) && b.isProperInterior(cands[i], tol)) {
            found = cands[i];
            ++nProper;
        }
    }
    if (nProper != 1) {
        return std::nullopt;
    }
    return Hit{found};
}

std::array<ExactCircularArc, 2>
ExactCircularArc::splitAt(const CoordinateXY& p) const
{
    const double tol = 1.0e-8 * std::max(1.0, m_r);
    if (!isProperInterior(p, tol)) {
        throw util::IllegalArgumentException(
            "splitAt requires a proper interior point of the arc");
    }
    const CoordinateXY mid0 = angularMid(m_start, p);
    const CoordinateXY mid1 = angularMid(p, m_end);
    return {
        ExactCircularArc(m_start, mid0, p),
        ExactCircularArc(p, mid1, m_end)
    };
}

std::optional<CookedArcs>
ExactCircularArc::cook(const ExactCircularArc& a, const ExactCircularArc& b)
{
    const auto h = hit(a, b);
    if (!h) {
        return std::nullopt;
    }
    const auto as = a.splitAt(h->point);
    const auto bs = b.splitAt(h->point);
    return CookedArcs{as[0], as[1], bs[0], bs[1], *h};
}

std::unique_ptr<CircularString>
ExactCircularArc::insertHit(const CircularString& cs,
                            std::size_t arcIndex,
                            const CoordinateXY& hitPt,
                            const CoordinateXY& mid0,
                            const CoordinateXY& mid1)
{
    const CoordinateSequence* pts = cs.getCoordinatesRO();
    const std::size_t i0 = arcIndex * 2;
    auto out = std::make_unique<CoordinateSequence>(
        pts->size() + 2, pts->hasZ(), pts->hasM());
    std::size_t dst = 0;
    auto copyPt = [&](std::size_t src) {
        CoordinateXYZM c;
        pts->getAt(src, c);
        out->setAt(c, dst++);
    };
    for (std::size_t i = 0; i <= i0; ++i) {
        copyPt(i);
    }
    out->setAt(mid0, dst++);
    out->setAt(hitPt, dst++);
    out->setAt(mid1, dst++);
    for (std::size_t i = i0 + 2; i < pts->size(); ++i) {
        copyPt(i);
    }
    return cs.getFactory()->createCircularString(std::move(out));
}

bool
ExactCircularArc::cook(const CircularString& a,
                       const CircularString& b,
                       std::unique_ptr<CircularString>& aOut,
                       std::unique_ptr<CircularString>& bOut,
                       Hit& hitOut)
{
    aOut.reset();
    bOut.reset();
    const auto& arcsA = a.getArcs();
    const auto& arcsB = b.getArcs();
    for (std::size_t i = 0; i < arcsA.size(); ++i) {
        ExactCircularArc ea(arcsA[i].p0(), arcsA[i].p1(), arcsA[i].p2());
        for (std::size_t j = 0; j < arcsB.size(); ++j) {
            ExactCircularArc eb(arcsB[j].p0(), arcsB[j].p1(), arcsB[j].p2());
            const auto cooked = cook(ea, eb);
            if (!cooked) {
                continue;
            }
            aOut = insertHit(a, i, cooked->hit.point,
                             cooked->a0.getMid(), cooked->a1.getMid());
            bOut = insertHit(b, j, cooked->hit.point,
                             cooked->b0.getMid(), cooked->b1.getMid());
            hitOut = cooked->hit;
            return true;
        }
    }
    return false;
}

} // namespace exactcurve
} // namespace algorithm
} // namespace geos

