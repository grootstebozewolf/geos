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

#include <geos/export.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Geometry.h>

#include <memory>

namespace geos {
namespace algorithm {
namespace exactcurve {

/**
 * Privileged ExactCurve primitive: one 3-control circular window.
 *
 * Maintainability: circumcircle lives here so callers do not copy the
 * determinant; colinear triples degrade to an exact chord.
 * Soundness: length() and pointAt() never densify.
 * Performance: static lengthOf is one hypot plus one multiply.
 * Port of JTS 9797c2c4.
 */
class GEOS_DLL ExactCircularArc {
public:
    ExactCircularArc(const geom::CoordinateXY& start,
                     const geom::CoordinateXY& mid,
                     const geom::CoordinateXY& end);

    static double lengthOf(const geom::CoordinateXY& start,
                           const geom::CoordinateXY& mid,
                           const geom::CoordinateXY& end);

    const geom::CoordinateXY& getStart() const { return m_start; }
    const geom::CoordinateXY& getMid() const { return m_mid; }
    const geom::CoordinateXY& getEnd() const { return m_end; }

    bool isArc() const { return m_arc; }
    bool isCcw() const { return m_ccw; }
    bool isExact() const { return true; }

    double radius() const { return m_r; }
    double sweep() const { return m_sweep; }
    double length() const;
    double chordLength() const { return m_start.distance(m_end); }

    geom::CoordinateXY center() const;

    geom::CoordinateXY pointAt(double t) const;
    std::unique_ptr<geom::Geometry> toLinear(double tolerance) const;

    bool chordLeArc() const;
    bool inArc(const geom::CoordinateXY& p, double radialTol) const;
    double circularSegmentArea() const;
    geom::CoordinateXY arcLengthCentroid() const;

    static bool tryCircumcircle(const geom::CoordinateXY& a,
                                const geom::CoordinateXY& b,
                                const geom::CoordinateXY& c,
                                double& cx, double& cy, double& r);

private:
    bool onSweep(const geom::CoordinateXY& p) const;
    static bool onSegment(const geom::CoordinateXY& p,
                          const geom::CoordinateXY& a,
                          const geom::CoordinateXY& b,
                          double tol);
    static int segmentCount(double radius, double sweep, double tolerance);
    static double directedSweep(double cx, double cy,
                                const geom::CoordinateXY& start,
                                const geom::CoordinateXY& mid,
                                const geom::CoordinateXY& end);
    static double ulp(double value);

    geom::CoordinateXY m_start;
    geom::CoordinateXY m_mid;
    geom::CoordinateXY m_end;
    double m_cx;
    double m_cy;
    double m_r;
    double m_a0;
    bool m_ccw;
    double m_sweep;
    bool m_arc;
};

} // namespace exactcurve
} // namespace algorithm
} // namespace geos
