/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2024-2026 ISciences, LLC
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/algorithm/Angle.h>
#include <geos/algorithm/CircularArcIntersector.h>
#include <geos/algorithm/CircularArcs.h>
#include <geos/algorithm/LineIntersector.h>
#include <geos/geom/CoordinateSequences.h>
#include <geos/math/DD.h>

#include <iomanip>
#include <tuple>

using geos::geom::CoordinateSequence;
using geos::geom::CoordinateXY;
using geos::geom::CircularArc;

namespace geos::algorithm {


// Interpolate the Z/M values of a point lying on the provided line segment
static void
interpolateSegmentZM(const CoordinateSequence& seq,
                     std::size_t ind0, std::size_t ind1,
                     CoordinateXY& pt, double& z, double& m)
{
    seq.applyAt(ind0, [&seq, &pt, ind1, &z, &m](const auto& p0) {
        using CoordinateType = std::decay_t<decltype(p0)>;

        const auto& p1 = seq.getAt<CoordinateType>(ind1);
        z = Interpolate::zGetOrInterpolate(pt, p0, p1);
        m = Interpolate::mGetOrInterpolate(pt, p0, p1);
    });
}


// Interpolate the Z/M values of an intersection point between two arcs
static void
interpolateZM(const CircularArc& arc0, const CircularArc& arc1, geom::CoordinateXYZM& pt)
{
    if (!std::isnan(pt.z) && !std::isnan(pt.m)) {
        return;
    }

    double z0, m0;
    double z1, m1;
    CircularArcs::interpolateZM(*arc0.getCoordinateSequence(), arc0.getCoordinatePosition(), arc0.getCenter(), arc0.isCCW(), pt, z0, m0);
    CircularArcs::interpolateZM(*arc1.getCoordinateSequence(), arc1.getCoordinatePosition(), arc1.getCenter(), arc1.isCCW(), pt, z1, m1);

    if (std::isnan(pt.z)) {
        pt.z = Interpolate::getOrAverage(z0, z1);
    }
    if (std::isnan(pt.m)) {
        pt.m = Interpolate::getOrAverage(m0, m1);
    }
}

// Interpolate the Z/M values of an intersection point between an arc and a segment
static void
interpolateZM(const CircularArc& arc0,
              const CoordinateSequence& seq,
              std::size_t ind0, std::size_t ind1,
              geom::CoordinateXYZM& pt)
{
    if (!std::isnan(pt.z) && !std::isnan(pt.m)) {
        return;
    }

    double z0, m0;
    double z1, m1;
    CircularArcs::interpolateZM(*arc0.getCoordinateSequence(), arc0.getCoordinatePosition(), arc0.getCenter(), arc0.isCCW(), pt, z0, m0);
    interpolateSegmentZM(seq, ind0, ind1, pt, z1, m1);

    if (std::isnan(pt.z)) {
        pt.z = Interpolate::getOrAverage(z0, z1);
    }
    if (std::isnan(pt.m)) {
        pt.m = Interpolate::getOrAverage(m0, m1);
    }
}


bool
CircularArcIntersector::hasIntersection(const geom::CoordinateXY &p) const {
    switch (nPt) {
        case 2: return intPt[1].equals2D(p) || intPt[0].equals2D(p);
        case 1: return intPt[0].equals2D(p);
        case 0: return false;
        default: break;
    }

    assert(0);
    return false;
}

static CoordinateXY&
closestPoint(CoordinateXY& p0, CoordinateXY& p1, int n, const CoordinateXY& q)
{
    if (n < 2) {
        return p0;
    }

    const double d0 = p0.distance(q);
    const double d1 = p1.distance(q);

    if (d0 < d1) {
        return p0;
    }

    return p1;
}

/// The circle through the arc's control points, in extended (double-double)
/// precision: a double-precision circumcenter carries a rounding error that
/// the intersection computation amplifies into the returned coordinates when
/// the circumradius is much larger than the coordinate magnitudes.
/// Orientation is canonicalized as in CircularArc::getCenter()/getRadius(),
/// so an arc and its reverse produce bitwise-identical circles.
static std::tuple<math::DD, math::DD, math::DD>
arcCircleDD(const CircularArc& arc)
{
    const CoordinateXY& q0 = arc.isCCW() ? arc.p0() : arc.p2();
    const CoordinateXY& q2 = arc.isCCW() ? arc.p2() : arc.p0();

    auto [cx, cy] = CircularArcs::getCenterDD(q0, arc.p1(), q2);

    math::DD rx = cx - q0.x;
    math::DD ry = cy - q0.y;

    return {cx, cy, rx*rx + ry*ry};
}

void
CircularArcIntersector::intersects(const CircularArc& arc, const CoordinateSequence& seq, std::size_t segPos0, std::size_t segPos1, bool useSegEndpoints)
{
    if (arc.isLinear()) {
        std::size_t arcPos0 = arc.getCoordinatePosition();
        intersects(*arc.getCoordinateSequence(), arcPos0, arcPos0 + 2,
                   seq, segPos0, segPos1);
        return;
    }

    reset();

    // TODO: envelope check?
    const CoordinateXY& ap0 = arc.p0<CoordinateXY>();
    const CoordinateXY& ap1 = arc.p1<CoordinateXY>();
    const CoordinateXY& ap2 = arc.p2<CoordinateXY>();
    const CoordinateXY& bp0 = seq.getAt<CoordinateXY>(segPos0);
    const CoordinateXY& bp1 = seq.getAt<CoordinateXY>(segPos1);

    // A segment endpoint equal to one of the arc's control points lies
    // exactly on both inputs and is always an intersection. Collect these
    // up front so that the floating-point circle-line computation below,
    // which can conclude "no intersection" for a grazing configuration,
    // cannot discard them.
    CoordinateXY shared[2];
    int nShared = 0;
    for (const CoordinateXY& a : {ap0, ap1, ap2}) {
        for (const CoordinateXY& b : {bp0, bp1}) {
            if (a.equals2D(b) && (nShared == 0 || !shared[0].equals2D(a))) {
                if (nShared < 2) {
                    shared[nShared++] = a;
                }
            }
        }
    }

    auto [cx, cy, rsq] = arcCircleDD(arc);

    CoordinateXY isect0, isect1;
    int nCand = CircularArcs::circleIntersectsLine(cx, cy, rsq, bp0, bp1, isect0, isect1);

    if (nCand == 2 && isect1.equals2D(isect0)) {
        nCand = 1;
    }

    // Check for exact endpoint-endpoint or endpoint-control point intersections
    // If found, replace the computed intersection points with an exact endpoint
    if (nCand > 0) {
        if (ap0 == bp0 || ap0 == bp1) {
            closestPoint(isect0, isect1, nCand, ap0) = ap0;
        }
        if (ap1 == bp0 || ap1 == bp1) {
            closestPoint(isect0, isect1, nCand, ap1) = ap1;
        }
        if (ap2 == bp0 || ap2 == bp1) {
            closestPoint(isect0, isect1, nCand, ap2) = ap2;
        }
    }

    Envelope segEnv(bp0, bp1);

    for (int i = 0; i < nShared; i++) {
        addArcSegmentIntersectionPoint(shared[i], arc, seq, segPos0, segPos1, useSegEndpoints);
    }

    if (nCand > 0 && nPt < 2 && !hasIntersection(isect0) && segEnv.contains(isect0) && arc.containsPointOnCircle(isect0)) {
        addArcSegmentIntersectionPoint(isect0, arc, seq, segPos0, segPos1, useSegEndpoints);
    }

    if (nCand > 1 && nPt < 2 && !hasIntersection(isect1) && segEnv.contains(isect1) && arc.containsPointOnCircle(isect1)) {
        addArcSegmentIntersectionPoint(isect1, arc, seq, segPos0, segPos1, useSegEndpoints);
    }

    switch (nPt) {
    case 2:
        result = TWO_POINT_INTERSECTION;
        break;
    case 1:
        result = ONE_POINT_INTERSECTION;
        break;
    default:
        result = NO_INTERSECTION;
    }
}

void
CircularArcIntersector::intersects(const CircularArc& arc1, const CircularArc& arc2)
{
    // Handle cases where one or both arcs are degenerate
    if (arc1.isLinear()) {
        if (arc2.isLinear()) {
            const auto arc1pos = arc1.getCoordinatePosition();
            const auto arc2pos = arc2.getCoordinatePosition();

            intersects(*arc1.getCoordinateSequence(), arc1pos, arc1pos + 2,
                       *arc2.getCoordinateSequence(), arc2pos, arc2pos + 2);
            return;
        } else {
            intersects(arc2, *arc1.getCoordinateSequence(), arc1.getCoordinatePosition(), arc1.getCoordinatePosition() + 2, true);
            return;
        }
    } else if (arc2.isLinear()) {
        intersects(arc1, *arc2.getCoordinateSequence(), arc2.getCoordinatePosition(), arc2.getCoordinatePosition() + 2, false);
        return;
    }

    reset();

    const CoordinateXY& ap0 = arc1.p0();
    const CoordinateXY& ap1 = arc1.p1();
    const CoordinateXY& ap2 = arc1.p2();
    const CoordinateXY& bp0 = arc2.p0();
    const CoordinateXY& bp1 = arc2.p1();
    const CoordinateXY& bp2 = arc2.p2();

    // A control point shared between the two arcs lies exactly on both of
    // them and is always an intersection. Collect these up front so that
    // the floating-point circle computation below, which can conclude
    // "no intersection" for circles that are tangent at the shared point,
    // cannot discard them.
    CoordinateXY shared[2];
    int nShared = 0;
    for (const CoordinateXY& a : {ap0, ap1, ap2}) {
        for (const CoordinateXY& b : {bp0, bp1, bp2}) {
            if (a.equals2D(b) && (nShared == 0 || !shared[0].equals2D(a))) {
                if (nShared < 2) {
                    shared[nShared++] = a;
                }
            }
        }
    }

    // Normalize arguments such that the computed intersection points do not depend
    // on the order of the input arcs
    const bool swapArgs = arc1.getCenter().compareTo(arc2.getCenter()) > 0;

    auto [c1x, c1y, r1sq] = arcCircleDD(swapArgs ? arc2 : arc1);
    auto [c2x, c2y, r2sq] = arcCircleDD(swapArgs ? arc1 : arc2);

    CoordinateXY isect0, isect1;
    int nCand = CircularArcs::circleIntersectsCircle(c1x, c1y, r1sq, c2x, c2y, r2sq, isect0, isect1);

    // TODO the coincidence test compares circles derived independently from each
    // arc's control points, so arcs on the same ideal circle that do not share
    // control points may still be reported as crossing rather than cocircular.
    if (nCand < 0) {
        computeCocircularIntersection(arc1, arc2);
    } else {
        if (nCand == 2 && isect1.equals2D(isect0)) {
            nCand = 1;
        }

        // Check to see if computed intersection points are inexact versions of an endpoint intersection
        if (nCand > 0) {
            if (ap0 == bp0 || ap0 == bp1 || ap0 == bp2) {
                closestPoint(isect0, isect1, nCand, ap0) = ap0;
            }
            if (ap1 == bp0 || ap1 == bp1 || ap1 == bp2) {
                closestPoint(isect0, isect1, nCand, ap1) = ap1;
            }
            if (ap2 == bp0 || ap2 == bp1 || ap2 == bp2) {
                closestPoint(isect0, isect1, nCand, ap2) = ap2;
            }
        }

        for (int i = 0; i < nShared; i++) {
            addArcArcIntersectionPoint(shared[i], arc1, arc2);
        }

        if (nCand > 0 && nPt < 2 && !hasIntersection(isect0) &&
            arc1.containsPointOnCircle(isect0) && arc2.containsPointOnCircle(isect0)) {
            addArcArcIntersectionPoint(isect0, arc1, arc2);
        }

        if (nCand > 1 && nPt < 2 && !hasIntersection(isect1) &&
            arc1.containsPointOnCircle(isect1) && arc2.containsPointOnCircle(isect1)) {
            addArcArcIntersectionPoint(isect1, arc1, arc2);
        }
    }

    if (nArc) {
        result = COCIRCULAR_INTERSECTION;
    }
    else {
        switch (nPt) {
        case 2:
            result = TWO_POINT_INTERSECTION;
            break;
        case 1:
            result = ONE_POINT_INTERSECTION;
            break;
        case 0:
            result = NO_INTERSECTION;
            break;
        default:
            assert(0);
        }
    }
}

void
CircularArcIntersector::intersects(const CoordinateSequence &p, std::size_t p0, std::size_t p1,
                                   const CoordinateSequence &q, std::size_t q0, std::size_t q1)
{
    LineIntersector li(precisionModel);
    li.computeIntersection(p, p0, p1, q, q0, q1);

    if (li.getIntersectionNum() == 2) {
        // FIXME this means a collinear intersection, so we should report as cocircular?
        intPt[0] = li.getIntersection(0);
        intPt[1] = li.getIntersection(1);
        result = TWO_POINT_INTERSECTION;
    } else if (li.getIntersectionNum() == 1) {
        intPt[0] = li.getIntersection(0);
        nPt = 1;
        result = ONE_POINT_INTERSECTION;
    } else {
        result = NO_INTERSECTION;
    }
}

/// Overwrite X/Y, and NaN Z/M values on the supplied point with those from the coordinate at the specified index
static void
setFromEndpoint(geom::CoordinateXYZM& pt, const CircularArc& arc, std::size_t index)
{
    arc.applyAt(index, [&pt](const auto& endpoint) {
        pt.x = endpoint.x;
        pt.y = endpoint.y;
        pt.z = Interpolate::zGet(pt, endpoint);
        pt.m = Interpolate::mGet(pt, endpoint);
    });
}

void
CircularArcIntersector::computeCocircularIntersection(const CircularArc& arc1, const CircularArc& arc2)
{
    const auto& center = arc1.getCenter();
    const double radius = arc1.getRadius();

    double ap0 = arc1.theta0();
    double ap1 = arc1.theta2();
    double bp0 = arc2.theta0();
    double bp1 = arc2.theta2();

    // Orientation of the result matches the first input
        bool resultArcIsCCW = true;

        // Make both inputs counter-clockwise for the purpose of determining intersections
        if (arc1.getOrientation() != Orientation::COUNTERCLOCKWISE) {
            std::swap(ap0, ap1);
            resultArcIsCCW = false;
        }
        if (arc2.getOrientation() != Orientation::COUNTERCLOCKWISE) {
            std::swap(bp0, bp1);
        }
        ap0 = Angle::normalizePositive(ap0);
        ap1 = Angle::normalizePositive(ap1);
        bp0 = Angle::normalizePositive(bp0);
        bp1 = Angle::normalizePositive(bp1);

        bool checkBp1inA = true;
        bool checkAcontained = true;

        // Possible intersection arrangements:
        // A contained within B
        // A overlaps B
        // B contained within A

        // check start of B within A?
        if (Angle::isWithinCCW(bp0, ap0, ap1)) {
            checkAcontained = false;
            const double start = bp0;
            const double end = Angle::nextCCW(start, bp1, ap1);

            if (end == bp1) {
                checkBp1inA = false;
            }

            if (start == end) {
                const CoordinateXY computedIntPt = CircularArcs::createPoint(center, radius, start);
                addArcArcIntersectionPoint(computedIntPt, arc1, arc2);
            }
            else {
                if (resultArcIsCCW) {
                    addCocircularIntersection(start, end, Orientation::COUNTERCLOCKWISE, arc1, arc2);
                }
                else {
                    addCocircularIntersection(end, start, Orientation::CLOCKWISE, arc1, arc2);
                }
            }
        }

        if (checkBp1inA && Angle::isWithinCCW(bp1, ap0, ap1)) {
            // end of B within A?
            checkAcontained = false;

            const double start = ap0;
            const double end = bp1;
            if (start == end) {
                const CoordinateXY computedIntPt = CircularArcs::createPoint(center, radius, start);
                addArcArcIntersectionPoint(computedIntPt, arc1, arc2);
            }
            else {
                if (resultArcIsCCW) {
                    addCocircularIntersection(start, end, Orientation::CLOCKWISE, arc1, arc2);
                }
                else {
                    addCocircularIntersection(end, start, Orientation::CLOCKWISE, arc1, arc2);
                }
            }
        }

        if (checkAcontained && Angle::isWithinCCW(ap0, bp0 , bp1) && ap0 != bp0  && ap0 != bp1 && Angle::isWithinCCW(ap1, bp0, bp1) && ap1 != bp1 && ap1 != bp0) {
            if (resultArcIsCCW) {
                addCocircularIntersection(ap0, ap1, Orientation::COUNTERCLOCKWISE, arc1, arc2);
            }
            else {
                addCocircularIntersection(ap1, ap0, Orientation::CLOCKWISE, arc1, arc2);
            }
        }
}

void
CircularArcIntersector::addCocircularIntersection(double startAngle, double endAngle, int orientation, const CircularArc& arc1, const CircularArc& arc2)
{
    const auto theta1  = CircularArcs::getMidpointAngle(startAngle, endAngle, orientation == Orientation::COUNTERCLOCKWISE);
    const CoordinateXY& center = arc1.getCenter();
    const double radius = arc1.getRadius();

    const bool constructZ = arc1.getCoordinateSequence()->hasZ() || arc2.getCoordinateSequence()->hasZ();
    const bool constructM = arc1.getCoordinateSequence()->hasM() || arc2.getCoordinateSequence()->hasM();

    CoordinateXYZM computedStartPt(CircularArcs::createPoint(center, radius, startAngle));
    CoordinateXYZM computedMidPt(CircularArcs::createPoint(center, radius, theta1));
    CoordinateXYZM computedEndPt(CircularArcs::createPoint(center, radius, endAngle));

    // Check to see if the endpoints of the intersection match the endpoints of either of
    // the input arcs. Use angles for the check to avoid missing an endpoint intersection from
    // inaccuracy in the point construction.
    if (startAngle == Angle::normalizePositive(arc1.theta0())) {
        computedStartPt = arc1.p0();
        setFromEndpoint(computedStartPt, arc1, 0);
    } else if (startAngle == Angle::normalizePositive(arc1.theta2())) {
        computedStartPt = arc1.p2();
        setFromEndpoint(computedStartPt, arc1, 2);
    } else if (startAngle == Angle::normalizePositive(arc2.theta0())) {
        computedStartPt = arc2.p0();
        setFromEndpoint(computedStartPt, arc2, 0);
    } else if (startAngle == Angle::normalizePositive(arc2.theta2())) {
        computedStartPt = arc2.p2();
        setFromEndpoint(computedStartPt, arc2, 2);
    }

    if (endAngle == Angle::normalizePositive(arc1.theta0())) {
        computedEndPt = arc1.p0();
        setFromEndpoint(computedEndPt, arc1, 0);
    } else if (endAngle == Angle::normalizePositive(arc1.theta2())) {
        computedEndPt = arc1.p2();
        setFromEndpoint(computedEndPt, arc1, 2);
    } else if (endAngle == Angle::normalizePositive(arc2.theta0())) {
        computedEndPt = arc2.p0();
        setFromEndpoint(computedEndPt, arc2, 0);
    } else if (endAngle == Angle::normalizePositive(arc2.theta2())) {
        computedEndPt = arc2.p2();
        setFromEndpoint(computedEndPt, arc2, 2);
    }

    interpolateZM(arc1, arc2, computedStartPt);
    interpolateZM(arc1, arc2, computedMidPt);
    interpolateZM(arc1, arc2, computedEndPt);

    if (precisionModel) {
        precisionModel->makePrecise(computedStartPt);
        precisionModel->makePrecise(computedMidPt);
        precisionModel->makePrecise(computedEndPt);
    }

    auto seq = std::make_unique<CoordinateSequence>(3, constructZ, constructM);
    seq->setAt(computedStartPt, 0);
    seq->setAt(computedMidPt, 1);
    seq->setAt(computedEndPt, 2);

    intArc[nArc++] = CircularArc(std::move(seq), 0, center, radius, orientation);
}

void
CircularArcIntersector::addArcArcIntersectionPoint(const CoordinateXY& computedIntPt, const CircularArc& arc1, const CircularArc& arc2) {
    CoordinateXYZM& newIntPt = intPt[nPt++];
    newIntPt = computedIntPt;

    if (precisionModel) {
        precisionModel->makePrecise(newIntPt);
    }

    if (computedIntPt.equals2D(arc1.p0())) {
        arc1.applyAt(0, [&newIntPt](const auto& endpoint) {
            newIntPt.z = Interpolate::zGet(newIntPt, endpoint);
            newIntPt.m = Interpolate::mGet(newIntPt, endpoint);
        });
    } else if (computedIntPt.equals2D(arc1.p2())) {
        arc1.applyAt(2, [&newIntPt](const auto& endpoint) {
            newIntPt.z = Interpolate::zGet(newIntPt, endpoint);
            newIntPt.m = Interpolate::mGet(newIntPt, endpoint);
        });
    }

    interpolateZM(arc1, arc2, newIntPt);
}

void
CircularArcIntersector::addArcSegmentIntersectionPoint(const CoordinateXY& computedIntPt, const CircularArc& arc1,
    const CoordinateSequence& seq, std::size_t pos0, std::size_t pos1, bool useSegEndpoints)
{
    CoordinateXYZM& newIntPt = intPt[nPt++];
    newIntPt = computedIntPt;

    if (precisionModel) {
        precisionModel->makePrecise(newIntPt);
    }

    for (int i = 0; i < 2; i++) {
        if (useSegEndpoints) {
            if (computedIntPt.equals2D(seq.getAt<CoordinateXY>(pos0))) {
                seq.applyAt(pos0, [&newIntPt](const auto& endpoint) {
                    newIntPt.z = Interpolate::zGet(newIntPt, endpoint);
                    newIntPt.m = Interpolate::mGet(newIntPt, endpoint);
                });
            }
            if (computedIntPt.equals2D(seq.getAt<CoordinateXY>(pos1))) {
                seq.applyAt(pos1, [&newIntPt](const auto& endpoint) {
                    newIntPt.z = Interpolate::zGet(newIntPt, endpoint);
                    newIntPt.m = Interpolate::mGet(newIntPt, endpoint);
                });
            }
        } else {
            if (computedIntPt.equals2D(arc1.p0())) {
                arc1.applyAt(0, [&newIntPt](const auto& endpoint) {
                    newIntPt.z = Interpolate::zGet(newIntPt, endpoint);
                    newIntPt.m = Interpolate::mGet(newIntPt, endpoint);
                });
            }
            if (computedIntPt.equals2D(arc1.p2())) {
                arc1.applyAt(2, [&newIntPt](const auto& endpoint) {
                    newIntPt.z = Interpolate::zGet(newIntPt, endpoint);
                    newIntPt.m = Interpolate::mGet(newIntPt, endpoint);
                });
            }
        }

        useSegEndpoints = !useSegEndpoints;
    }

    interpolateZM(arc1, seq, pos0, pos1, newIntPt);
}


}