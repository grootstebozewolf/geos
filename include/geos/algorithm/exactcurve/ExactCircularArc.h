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

#include <geos/algorithm/exactcurve/ExactCurve.h>
#include <geos/export.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Geometry.h>

#include <array>
#include <memory>
#include <optional>

namespace geos {
namespace geom {
class CircularString;
}
namespace algorithm {
namespace exactcurve {

/**
 * Unique proper intersection of two circular arcs.
 *
 * Empty when the arcs miss, touch only at an endpoint, are cocircular,
 * or meet at more than one interior point.
 */
struct Hit {
    geom::CoordinateXY point;
};

class ExactCircularArc;
struct CookedArcs;

/**
 * Closed-form circular arc through three control points.
 *
 * Collinear triples degrade to a chord. length() and pointAt() use the
 * circumcircle (radius times sweep) and do not densify. toLinear() is
 * the only path that emits a LineString.
 *
 * This type is XY-only: Z and M are not stored or interpolated.
 * cook(CircularString) copies existing vertices with their source
 * ordinates; the inserted hit and angular midpoints are XY.
 *
 * Inspired by JTS CircularArc (commit 9797c2c4).
 */
class GEOS_DLL ExactCircularArc : public ExactCurve {
public:
    ExactCircularArc(const geom::CoordinateXY& start,
                     const geom::CoordinateXY& mid,
                     const geom::CoordinateXY& end);

    static double lengthOf(const geom::CoordinateXY& start,
                           const geom::CoordinateXY& mid,
                           const geom::CoordinateXY& end);

    const geom::CoordinateXY& getStart() const override { return m_start; }
    const geom::CoordinateXY& getMid() const { return m_mid; }
    const geom::CoordinateXY& getEnd() const override { return m_end; }

    bool isArc() const { return m_arc; }
    bool isCcw() const { return m_ccw; }
    bool isExact() const override { return true; }

    double radius() const { return m_r; }
    double sweep() const { return m_sweep; }
    double length() const override;
    double chordLength() const { return m_start.distance(m_end); }

    geom::CoordinateXY center() const;

    geom::CoordinateXY pointAt(double t) const override;
    std::unique_ptr<geom::Geometry> toLinear(double tolerance) const override;

    bool chordLeArc() const;
    bool inArc(const geom::CoordinateXY& p, double radialTol) const;
    double circularSegmentArea() const;
    geom::CoordinateXY arcLengthCentroid() const;

    static bool tryCircumcircle(const geom::CoordinateXY& a,
                                const geom::CoordinateXY& b,
                                const geom::CoordinateXY& c,
                                double& cx, double& cy, double& r);

    /**
     * Unique proper circular-circular intersection, or none.
     * Does not call getLinearized / addLinearizedPoints / Densifier.
     */
    static std::optional<Hit> hit(const ExactCircularArc& a,
                                  const ExactCircularArc& b);

    /**
     * Insert an interior vertex, returning two child windows.
     * Midpoints are angular midpoints on this circumcircle.
     */
    std::array<ExactCircularArc, 2> splitAt(const geom::CoordinateXY& p) const;

    /**
     * If the arcs have a unique proper intersection (a hit), split each
     * at that point and return the four children plus the hit.
     * isExact() remains true on every child. Does not densify.
     */
    static std::optional<CookedArcs> cook(const ExactCircularArc& a,
                                          const ExactCircularArc& b);

    /**
     * Same as cook(arc, arc) for the first pair of 3-point windows that
     * have a unique proper intersection. Each output CircularString is
     * the input with one vertex inserted (3 controls become 5).
     */
    static bool cook(const geom::CircularString& a,
                     const geom::CircularString& b,
                     std::unique_ptr<geom::CircularString>& aOut,
                     std::unique_ptr<geom::CircularString>& bOut,
                     Hit& hit);

private:
    bool onSweep(const geom::CoordinateXY& p) const;
    bool isProperInterior(const geom::CoordinateXY& p, double radialTol) const;
    geom::CoordinateXY angularMid(const geom::CoordinateXY& from,
                                  const geom::CoordinateXY& to) const;
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
    static std::unique_ptr<geom::CircularString> insertHit(
        const geom::CircularString& cs,
        std::size_t arcIndex,
        const geom::CoordinateXY& hitPt,
        const geom::CoordinateXY& mid0,
        const geom::CoordinateXY& mid1);

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

/**
 * Four child arcs after inserting a hit vertex on each input.
 */
struct CookedArcs {
    ExactCircularArc a0;
    ExactCircularArc a1;
    ExactCircularArc b0;
    ExactCircularArc b1;
    Hit hit;
};

} // namespace exactcurve
} // namespace algorithm
} // namespace geos
