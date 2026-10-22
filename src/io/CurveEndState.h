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

#pragma once

#include <geos/constants.h>
#include <geos/geom/CircularString.h>
#include <geos/geom/Clothoid.h>

#include <cmath>

namespace geos {
namespace io {

/// Largest difference between a short-form CLOTHOID k0 and the previous
/// member's end curvature. It covers the ProRail arc-radius rounding and
/// rejects k0 = 0.02 after a line.
constexpr double CLOTHOID_START_CURVATURE_TOLERANCE = 1e-7;

/// Heading and curvature at the end of a simple curve. Positive curvature is CCW.
inline bool
curveEndState(const geom::SimpleCurve& curve, double& heading, double& curvature)
{
    if (curve.isEmpty() || curve.getNumPoints() < 2) {
        return false;
    }

    if (curve.getGeometryTypeId() == geom::GEOS_CLOTHOID) {
        const auto& clothoid = static_cast<const geom::Clothoid&>(curve);
        heading = clothoid.getEndHeading();
        curvature = clothoid.getEndCurvature();
        return true;
    }

    if (curve.getGeometryTypeId() == geom::GEOS_CIRCULARSTRING) {
        const auto& circular = static_cast<const geom::CircularString&>(curve);
        if (circular.getNumPoints() < 3) {
            return false;
        }
        const auto& arcs = circular.getArcs();
        if (arcs.empty()) {
            return false;
        }
        const geom::CircularArc& arc = arcs.back();
        if (!arc.isLinear() && std::isfinite(arc.getRadius()) && arc.getRadius() > 0.0) {
            const geom::CoordinateXY& center = arc.getCenter();
            const geom::CoordinateXY& end = arc.p2();
            const double rx = end.x - center.x;
            const double ry = end.y - center.y;
            if (arc.isCCW()) {
                heading = std::atan2(rx, -ry);
                curvature = 1.0 / arc.getRadius();
            } else {
                heading = std::atan2(-rx, ry);
                curvature = -1.0 / arc.getRadius();
            }
            return std::isfinite(heading) && std::isfinite(curvature);
        }
    }

    const geom::CoordinateSequence* seq = curve.getCoordinatesRO();
    const geom::CoordinateXY& a = seq->getAt<geom::CoordinateXY>(seq->size() - 2);
    const geom::CoordinateXY& b = seq->getAt<geom::CoordinateXY>(seq->size() - 1);
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    if (!(std::hypot(dx, dy) > 0.0)) {
        return false;
    }
    heading = std::atan2(dy, dx);
    curvature = 0.0;
    return true;
}

/// Short-form WKT is safe only when parsing it again rebuilds this clothoid.
inline bool
clothoidContinues(const geom::SimpleCurve& previous, const geom::Clothoid& clothoid)
{
    if (clothoid.isEmpty() || previous.isEmpty() || !(clothoid.getScaleFactor() > 0.0)) {
        return false;
    }
    if (clothoid.getStartCoordinate().distance(previous.getEndCoordinate()) > 1e-8) {
        return false;
    }
    const geom::CoordinateSequence* clothoidSeq = clothoid.getCoordinatesRO();
    const geom::CoordinateSequence* previousSeq = previous.getCoordinatesRO();
    if (clothoid.hasZ()) {
        if (clothoidSeq->getZ(0) != clothoidSeq->getZ(clothoidSeq->size() - 1)) {
            return false;
        }
    }
    if (clothoid.hasZ() || previous.hasZ()) {
        const double clothoidZ = clothoid.hasZ() ? clothoidSeq->getZ(0) : geos::DoubleNotANumber;
        const double previousZ = previous.hasZ()
            ? previousSeq->getZ(previousSeq->size() - 1) : geos::DoubleNotANumber;
        if (!(clothoidZ == previousZ)) {
            return false;
        }
    }
    if (clothoid.hasM() && clothoid.getStartM() != clothoid.getEndM()) {
        return false;
    }
    if (clothoid.hasM() || previous.hasM()) {
        const double clothoidM = clothoid.hasM() ? clothoidSeq->getM(0) : geos::DoubleNotANumber;
        const double previousM = previous.hasM()
            ? previousSeq->getM(previousSeq->size() - 1) : geos::DoubleNotANumber;
        if (!(clothoidM == previousM)) {
            return false;
        }
    }

    const geom::ClothoidReferenceLocation loc = clothoid.getReferenceLocation();
    const double u = std::hypot(loc.u.x, loc.u.y, loc.u.z);
    const double v = std::hypot(loc.v.x, loc.v.y, loc.v.z);
    const double dot = loc.u.x * loc.v.x + loc.u.y * loc.v.y + loc.u.z * loc.v.z;
    if (std::abs(u - 1.0) > 1e-8 || std::abs(v - 1.0) > 1e-8 || std::abs(dot) > 1e-8) {
        return false;
    }

    double heading = 0.0;
    double curvature = 0.0;
    if (!curveEndState(previous, heading, curvature)) {
        return false;
    }
    const double turn = std::remainder(clothoid.getStartHeading() - heading, 2.0 * MATH_PI);
    if (!(std::abs(turn) <= 1e-8)) {
        return false;
    }
    return std::abs(clothoid.getStartCurvature() - curvature) <= CLOTHOID_START_CURVATURE_TOLERANCE;
}

}
}
