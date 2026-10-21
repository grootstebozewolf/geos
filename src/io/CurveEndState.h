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

/// Heading and curvature at the end of a simple curve. Positive curvature is CCW.
inline bool
curveEndState(const geom::SimpleCurve& curve, double& heading, double& curvature)
{
    if (curve.isEmpty() || curve.getNumPoints() < 2) {
        return false;
    }

    if (curve.getGeometryTypeId() == geom::GEOS_CLOTHOID) {
        const auto& clothoid = static_cast<const geom::Clothoid&>(curve);
        heading = clothoid.getEndTangent();
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
    if (clothoid.hasZ()) {
        const geom::CoordinateSequence* seq = clothoid.getCoordinatesRO();
        if (seq->getZ(0) != seq->getZ(seq->size() - 1)) {
            return false;
        }
    }
    if (clothoid.hasM() && clothoid.getStartM() != clothoid.getEndM()) {
        return false;
    }

    const geom::ClothoidReferenceLocation loc = clothoid.getReferenceLocation();
    const double u = std::hypot(loc.ux, loc.uy, loc.uz);
    const double v = std::hypot(loc.vx, loc.vy, loc.vz);
    const double dot = loc.ux * loc.vx + loc.uy * loc.vy + loc.uz * loc.vz;
    if (std::abs(u - 1.0) > 1e-8 || std::abs(v - 1.0) > 1e-8 || std::abs(dot) > 1e-8) {
        return false;
    }

    double heading = 0.0;
    double curvature = 0.0;
    if (!curveEndState(previous, heading, curvature)) {
        return false;
    }
    const double turn = std::remainder(clothoid.getStartTangent() - heading, 2.0 * MATH_PI);
    if (!(std::abs(turn) <= 1e-8)) {
        return false;
    }
    const double scale = std::max(1.0, std::max(std::abs(curvature), std::abs(clothoid.getStartCurvature())));
    return std::abs(clothoid.getStartCurvature() - curvature) <= 1e-6 * scale;
}

}
}
