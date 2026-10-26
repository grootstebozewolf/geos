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

#include <geos/geom/Coordinate.h>
#include <geos/geom/Envelope.h>

#include <cstddef>
#include <vector>

namespace geos {
namespace geom {
class CircularArc;
class Clothoid;
}
namespace algorithm {
namespace clothoid {

/// The plan (XY) image of a non-empty clothoid, evaluated by constructive
/// distance t between its start and end distances.
class PlacedClothoid {
public:
    explicit PlacedClothoid(const geom::Clothoid& curve);

    double startDistance() const
    {
        return m_t0;
    }

    double endDistance() const
    {
        return m_t1;
    }

    /// Point at constructive distance t. The ends are the stored endpoints.
    geom::CoordinateXY pointAt(double t) const;

    /// dP/dt at constructive distance t.
    geom::CoordinateXY tangentAt(double t) const;

    /// dT/dtheta at constructive distance t.
    geom::CoordinateXY normalAt(double t) const;

    /// d theta / dt at constructive distance t.
    double headingRate(double t) const;

    /// A point within this distance of the curve is on it: the endpoint snap
    /// tolerance, 32 ulps of the reference location's extent.
    double tolerance() const
    {
        return m_tolerance;
    }

    double maxSingularValue() const
    {
        return m_sigmaMax;
    }

    const geom::Envelope& envelope() const
    {
        return m_envelope;
    }

    /// Constructive distances in travel order, from start to end, between
    /// which x and y are both monotone, so the box of a span's endpoints holds it.
    const std::vector<double>& monotoneBreaks() const
    {
        return m_breaks;
    }

    /// Constructive distances strictly between the start and end distances,
    /// ascending, where the tangent is perpendicular to (dx, dy).
    std::vector<double> perpendicularTangents(double dx, double dy) const;

private:
    double m_ox;
    double m_oy;
    double m_ux;
    double m_uy;
    double m_vx;
    double m_vy;
    double m_scale;
    double m_t0;
    double m_t1;
    geom::CoordinateXY m_start;
    geom::CoordinateXY m_end;
    double m_sigmaMax;
    double m_tolerance;
    geom::Envelope m_envelope;
    std::vector<double> m_breaks;
};

/// Closest points between a clothoid and another primitive.
struct Nearest {
    double distance;
    geom::CoordinateXY onCurve;
    geom::CoordinateXY onTarget;
};

Nearest nearestToPoint(const PlacedClothoid& curve, const geom::CoordinateXY& q);

Nearest nearestToSegment(const PlacedClothoid& curve,
                         const geom::CoordinateXY& p0, const geom::CoordinateXY& p1);

Nearest nearestToArc(const PlacedClothoid& curve, const geom::CircularArc& arc);

Nearest nearestToClothoid(const PlacedClothoid& curve, const PlacedClothoid& other);

/// Counts crossings of the ray from q towards +x with the half-open vertex
/// rule of RayCrossingCounter::countSegment. Sets onCurve when q is within
/// tolerance of the curve; crossings are then not counted.
void countRayCrossings(const PlacedClothoid& curve, const geom::CoordinateXY& q,
                       std::size_t& crossings, bool& onCurve);

}
}
}
