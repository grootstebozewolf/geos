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
#include <geos/algorithm/Orientation.h>
#include <geos/constants.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/util/IllegalArgumentException.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

using geos::geom::CoordinateSequence;
using geos::geom::CoordinateXY;
using geos::geom::Geometry;
using geos::geom::GeometryFactory;

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
    std::uint64_t bits = 0;
    std::memcpy(&bits, &abs, sizeof(bits));
    bits += 1;
    double next;
    std::memcpy(&next, &bits, sizeof(next));
    return next - abs;
}

} // namespace exactcurve
} // namespace algorithm
} // namespace geos
