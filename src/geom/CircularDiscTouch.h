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

#pragma once

#include <geos/algorithm/Orientation.h>
#include <geos/constants.h>
#include <geos/geom/CircularString.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Curve.h>
#include <geos/geom/CurvePolygon.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/IntersectionMatrix.h>

#include <array>
#include <cmath>
#include <memory>
#include <optional>

namespace geos {
namespace geom {
namespace disc_touch {

using geos::algorithm::Orientation;

constexpr double TWO_PI = 2.0 * MATH_PI;
constexpr double SWEEP_EPS = 1.0e-9;
constexpr const char* AREA_EXT_TANGENT = "FF2F01212";

inline double
normPos(double angle)
{
    angle = std::fmod(angle, TWO_PI);
    if (angle < 0.0) {
        angle += TWO_PI;
    }
    return angle;
}

inline std::optional<std::array<double, 3>>
circumcircle(const CoordinateXY& a, const CoordinateXY& b, const CoordinateXY& c)
{
    if (Orientation::index(a, b, c) == Orientation::COLLINEAR) {
        return std::nullopt;
    }
    const double d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if (d == 0.0) {
        return std::nullopt;
    }
    const double a2 = a.x * a.x + a.y * a.y;
    const double b2 = b.x * b.x + b.y * b.y;
    const double c2 = c.x * c.x + c.y * c.y;
    const double ux = (a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / d;
    const double uy = (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / d;
    const double r = std::hypot(a.x - ux, a.y - uy);
    if (!std::isfinite(r) || r == 0.0) {
        return std::nullopt;
    }
    return std::array<double, 3>{ux, uy, r};
}

inline double
signedSweep(const CoordinateXY& start, const CoordinateXY& mid,
            const CoordinateXY& end, const std::array<double, 3>& c)
{
    const double a0 = std::atan2(start.y - c[1], start.x - c[0]);
    const double aMid = std::atan2(mid.y - c[1], mid.x - c[0]);
    const double a1 = std::atan2(end.y - c[1], end.x - c[0]);
    const bool ccw = normPos(aMid - a0) < normPos(a1 - a0);
    double sweep = ccw ? normPos(a1 - a0) : -normPos(a0 - a1);
    if (sweep == 0.0) {
        sweep = ccw ? TWO_PI : -TWO_PI;
    }
    return sweep;
}

inline std::optional<std::array<double, 3>>
fullCircle(const Geometry& ring)
{
    if (ring.getGeometryTypeId() != GEOS_CIRCULARSTRING
            || ring.isEmpty() || ring.getNumPoints() < 5) {
        return std::nullopt;
    }
    const auto* cs = dynamic_cast<const CircularString*>(&ring);
    if (!cs || !cs->isClosed()) {
        return std::nullopt;
    }
    std::unique_ptr<CoordinateSequence> seq(cs->getCoordinates());
    std::optional<std::array<double, 3>> found;
    double sweep = 0.0;
    for (std::size_t i = 0; i + 2 < seq->size(); i += 2) {
        auto c = circumcircle(seq->getAt<CoordinateXY>(i),
                              seq->getAt<CoordinateXY>(i + 1),
                              seq->getAt<CoordinateXY>(i + 2));
        if (!c) {
            return std::nullopt;
        }
        if (!found) {
            found = c;
        }
        else if (std::hypot((*found)[0] - (*c)[0], (*found)[1] - (*c)[1]) > 1.0e-9
                || std::fabs((*found)[2] - (*c)[2]) > 1.0e-9) {
            return std::nullopt;
        }
        sweep += signedSweep(seq->getAt<CoordinateXY>(i),
                             seq->getAt<CoordinateXY>(i + 1),
                             seq->getAt<CoordinateXY>(i + 2), *c);
    }
    if (!found || std::fabs(std::fabs(sweep) - TWO_PI) > SWEEP_EPS) {
        return std::nullopt;
    }
    return found;
}

inline std::optional<std::array<double, 3>>
circularDisc(const Geometry& g)
{
    if (g.getGeometryTypeId() != GEOS_CURVEPOLYGON || g.isEmpty()) {
        return std::nullopt;
    }
    const auto* cp = dynamic_cast<const CurvePolygon*>(&g);
    if (!cp || cp->getNumInteriorRing() > 0) {
        return std::nullopt;
    }
    const Curve* shell = cp->getExteriorRing();
    if (shell == nullptr) {
        return std::nullopt;
    }
    return fullCircle(*shell);
}

inline std::unique_ptr<IntersectionMatrix>
tryDiscExternalTouch(const Geometry& a, const Geometry& b)
{
    auto da = circularDisc(a);
    auto db = circularDisc(b);
    if (!da || !db) {
        return nullptr;
    }
    const double dx = (*da)[0] - (*db)[0];
    const double dy = (*da)[1] - (*db)[1];
    const double d2 = dx * dx + dy * dy;
    const double sum = (*da)[2] + (*db)[2];
    const double sum2 = sum * sum;
    if (d2 != sum2) {
        return nullptr;
    }
    return std::make_unique<IntersectionMatrix>(AREA_EXT_TANGENT);
}

} // namespace disc_touch
} // namespace geom
} // namespace geos
