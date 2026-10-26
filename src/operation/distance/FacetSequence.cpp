/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2016 Daniel Baston
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: operation/distance/FacetSequence.java (aff11591)
 *
 **********************************************************************/

#include <algorithm/ClothoidGeometry.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/LineSegment.h>
#include <geos/algorithm/CircularArcs.h>
#include <geos/algorithm/Distance.h>
#include <geos/geom/CircularArc.h>
#include <geos/geom/Clothoid.h>
#include <geos/operation/distance/FacetSequence.h>

#include <memory>

using namespace geos::geom;
using namespace geos::operation::distance;
using namespace geos::algorithm;


FacetSequence::FacetSequence(const CoordinateSequence* p_pts, std::size_t p_start, std::size_t p_end, bool p_arcs)
    : pts(p_pts)
    , start(p_start)
    , end(p_end)
    , arcs(p_arcs)
    , clothoid(nullptr)
{
    computeEnvelope();
}

FacetSequence::FacetSequence(const Clothoid* p_clothoid)
    : pts(p_clothoid->getCoordinatesRO())
    , start(0)
    , end(p_clothoid->getCoordinatesRO()->size())
    , arcs(false)
    , clothoid(p_clothoid)
{
    computeEnvelope();
}

size_t
FacetSequence::size() const
{
    return end - start;
}

bool
FacetSequence::isPoint() const
{
    return end - start == 1;
}

double
FacetSequence::distance(const FacetSequence& facetSeq) const
{
    if (arcs || facetSeq.arcs || clothoid || facetSeq.clothoid) {
        return computeDistanceWithArcs(facetSeq, nullptr);
    }

    bool isPointThis = isPoint();
    bool isPointOther = facetSeq.isPoint();

    if (isPointThis && isPointOther) {
        const Coordinate& pt = pts->getAt(start);
        const Coordinate& seqPt = facetSeq.pts->getAt(facetSeq.start);
        return pt.distance(seqPt);
    }
    else if (isPointThis) {
        const Coordinate& pt = pts->getAt(start);
        return computeDistancePointLine(pt, facetSeq, nullptr);
    }
    else if (isPointOther) {
        const Coordinate& seqPt = facetSeq.pts->getAt(facetSeq.start);
        return computeDistancePointLine(seqPt, *this, nullptr);
    }
    else {
        return computeDistanceLineLine(facetSeq, nullptr);
    }
}

/*
* Rather than get bent out of shape about returning a pointer
* just return the whole mess, since it only ends up holding two
* locations.
*/
std::vector<Coordinate>
FacetSequence::nearestLocations(const FacetSequence& facetSeq)  const
{
    std::vector<Coordinate> locs;
    if (arcs || facetSeq.arcs || clothoid || facetSeq.clothoid) {
        computeDistanceWithArcs(facetSeq, &locs);
        return locs;
    }

    bool isPointThis = isPoint();
    bool isPointOther = facetSeq.isPoint();
    if (isPointThis && isPointOther) {
        const Coordinate& pt = pts->getAt(start);
        const Coordinate& seqPt = facetSeq.pts->getAt(facetSeq.start);
        locs.clear();
        locs.push_back(pt);
        locs.push_back(seqPt);
    }
    else if (isPointThis) {
        const Coordinate& pt = pts->getAt(start);
        computeDistancePointLine(pt, facetSeq, &locs);
    }
    else if (isPointOther) {
        const Coordinate& seqPt = facetSeq.pts->getAt(facetSeq.start);
        computeDistancePointLine(seqPt, *this, &locs);
        // unflip the locations
        Coordinate tmp = locs[0];
        locs[0] = locs[1];
        locs[1] = tmp;
    }
    else {
        computeDistanceLineLine(facetSeq, &locs);
    }
    return locs;
}

double
FacetSequence::computeDistancePointLine(const Coordinate& pt,
                                        const FacetSequence& facetSeq,
                                        std::vector<Coordinate> *locs) const
{
    double minDistance = DoubleInfinity;

    for (std::size_t i = facetSeq.start; i < facetSeq.end - 1; i++) {
        const Coordinate& q0 = facetSeq.pts->getAt(i);
        const Coordinate& q1 = facetSeq.pts->getAt(i + 1);
        double dist = Distance::pointToSegment(pt, q0, q1);
        if (dist < minDistance || (locs != nullptr && locs->empty())) {
            minDistance = dist;
            if (locs != nullptr) {
                updateNearestLocationsPointLine(pt, q0, q1, locs);
            }
            if (minDistance <= 0.0) {
                return minDistance;
            }
        }
    }

    return minDistance;
}

void
FacetSequence::updateNearestLocationsPointLine(
    const Coordinate& pt, const Coordinate& q0, const Coordinate &q1,
    std::vector<Coordinate> *locs) const
{
    LineSegment seg(q0, q1);
    Coordinate segClosestPoint;
    seg.closestPoint(pt, segClosestPoint);
    locs->clear();
    locs->push_back(pt);
    locs->push_back(segClosestPoint);
    return;
}

double
FacetSequence::computeDistanceLineLine(const FacetSequence& facetSeq, std::vector<Coordinate> *locs) const
{
    double minDistance = DoubleInfinity;

    for(std::size_t i = start; i < end - 1; i++) {
        const Coordinate& p0 = pts->getAt(i);
        const Coordinate& p1 = pts->getAt(i + 1);

        // Avoid calculating distance from zero-length segment
        if (p0 == p1)
            continue;

        Envelope pEnv(p0, p1);
        if (pEnv.distanceSquared(*facetSeq.getEnvelope()) > minDistance*minDistance) {
            continue;
        }

        for(std::size_t j = facetSeq.start; j < facetSeq.end - 1; j++) {
            const Coordinate& q0 = facetSeq.pts->getAt(j);
            const Coordinate& q1 = facetSeq.pts->getAt(j + 1);

            // Avoid calculating distance to zero-length segment
            if (q0 == q1)
                continue;

            Envelope qEnv(q0, q1);
            if (pEnv.distanceSquared(qEnv) > minDistance*minDistance) {
                continue;
            }

            double dist = Distance::segmentToSegment(p0, p1, q0, q1);
            if(dist <= minDistance) {
                minDistance = dist;
                if(locs != nullptr) {
                    updateNearestLocationsLineLine(p0, p1, q0, q1, locs);
                }
                if(minDistance <= 0.0) return minDistance;
            }
        }
    }

    return minDistance;
}

void
FacetSequence::updateNearestLocationsLineLine(
    const Coordinate& p0, const Coordinate& p1,
    const Coordinate& q0, const Coordinate &q1,
    std::vector<Coordinate> *locs) const
{
    LineSegment seg0(p0, p1);
    LineSegment seg1(q0, q1);

    auto closestPts = seg0.closestPoints(seg1);

    locs->clear();
    locs->push_back(closestPts[0]);
    locs->push_back(closestPts[1]);
}

/*
* The facets of a sequence are its point (if it has a single coordinate),
* its segments, or its arcs.
*/
std::size_t
FacetSequence::numElements() const
{
    if (clothoid || isPoint()) {
        return 1;
    }
    return arcs ? (end - start - 1) / 2 : end - start - 1;
}

Envelope
FacetSequence::elementEnvelope(std::size_t i) const
{
    Envelope e;
    if (clothoid) {
        e = env;
    }
    else if (isPoint()) {
        e.expandToInclude(pts->getAt<CoordinateXY>(start));
    }
    else if (arcs) {
        const std::size_t k = start + 2*i;
        CircularArcs::expandEnvelope(e, pts->getAt<CoordinateXY>(k),
                                     pts->getAt<CoordinateXY>(k + 1),
                                     pts->getAt<CoordinateXY>(k + 2));
    }
    else {
        e.expandToInclude(pts->getAt<CoordinateXY>(start + i));
        e.expandToInclude(pts->getAt<CoordinateXY>(start + i + 1));
    }
    return e;
}

std::array<CoordinateXY, 2>
FacetSequence::closestPoints(std::size_t i, const FacetSequence& facetSeq, std::size_t j) const
{
    if (facetSeq.clothoid && !clothoid) {
        auto pts2 = facetSeq.closestPoints(j, *this, i);
        return {pts2[1], pts2[0]};
    }

    if (clothoid) {
        using namespace geos::algorithm::clothoid;
        const PlacedClothoid placed(*clothoid);
        Nearest n;
        if (facetSeq.clothoid) {
            n = nearestToClothoid(placed, PlacedClothoid(*facetSeq.clothoid));
        }
        else if (facetSeq.isPoint()) {
            n = nearestToPoint(placed, facetSeq.pts->getAt<CoordinateXY>(facetSeq.start));
        }
        else if (facetSeq.arcs) {
            n = nearestToArc(placed, CircularArc(*facetSeq.pts, facetSeq.start + 2*j));
        }
        else {
            n = nearestToSegment(placed, facetSeq.pts->getAt<CoordinateXY>(facetSeq.start + j),
                                 facetSeq.pts->getAt<CoordinateXY>(facetSeq.start + j + 1));
        }
        return {n.onCurve, n.onTarget};
    }

    if (isPoint()) {
        const CoordinateXY& p = pts->getAt<CoordinateXY>(start);
        if (facetSeq.isPoint()) {
            return {p, facetSeq.pts->getAt<CoordinateXY>(facetSeq.start)};
        }
        if (facetSeq.arcs) {
            return {p, CircularArc(*facetSeq.pts, facetSeq.start + 2*j).closestPoint(p)};
        }
        return {p, LineSegment::closestPoint(facetSeq.pts->getAt<CoordinateXY>(facetSeq.start + j),
                                             facetSeq.pts->getAt<CoordinateXY>(facetSeq.start + j + 1), p)};
    }

    if (facetSeq.isPoint() || !arcs) {
        auto pts2 = facetSeq.closestPoints(j, *this, i);
        return {pts2[1], pts2[0]};
    }

    const CircularArc arc(*pts, start + 2*i);
    if (facetSeq.arcs) {
        return arc.closestPoints(CircularArc(*facetSeq.pts, facetSeq.start + 2*j));
    }
    return arc.closestPoints(facetSeq.pts->getAt<CoordinateXY>(facetSeq.start + j),
                             facetSeq.pts->getAt<CoordinateXY>(facetSeq.start + j + 1));
}

double
FacetSequence::computeDistanceWithArcs(const FacetSequence& facetSeq, std::vector<Coordinate>* locs) const
{
    double minDistance = DoubleInfinity;

    for (std::size_t i = 0; i < numElements(); i++) {
        const Envelope pEnv = elementEnvelope(i);
        if (pEnv.distanceSquared(*facetSeq.getEnvelope()) > minDistance*minDistance) {
            continue;
        }

        for (std::size_t j = 0; j < facetSeq.numElements(); j++) {
            if (pEnv.distanceSquared(facetSeq.elementEnvelope(j)) > minDistance*minDistance) {
                continue;
            }

            const auto closestPts = closestPoints(i, facetSeq, j);
            const double dist = closestPts[0].distance(closestPts[1]);
            if (dist < minDistance || (locs != nullptr && locs->empty())) {
                minDistance = dist;
                if (locs != nullptr) {
                    locs->clear();
                    locs->emplace_back(closestPts[0]);
                    locs->emplace_back(closestPts[1]);
                }
                if (minDistance <= 0.0) {
                    return minDistance;
                }
            }
        }
    }

    return minDistance;
}

void
FacetSequence::computeEnvelope()
{
    env = Envelope();
    if (clothoid) {
        // A point within the evaluation tolerance of the clothoid is on it.
        env = *clothoid->getEnvelopeInternal();
        env.expandBy(algorithm::clothoid::PlacedClothoid(*clothoid).tolerance());
        return;
    }
    if (arcs) {
        for (std::size_t i = start; i + 2 < end; i += 2) {
            CircularArcs::expandEnvelope(env, pts->getAt<CoordinateXY>(i),
                                         pts->getAt<CoordinateXY>(i + 1),
                                         pts->getAt<CoordinateXY>(i + 2));
        }
        return;
    }
    for(std::size_t i = start; i < end; i++) {
        env.expandToInclude(pts->getAt(i));
    }
}

const Envelope*
FacetSequence::getEnvelope() const
{
    return &env;
}

const Coordinate*
FacetSequence::getCoordinate(std::size_t index) const
{
    return &(pts->getAt(start + index));
}

std::size_t
FacetSequence::normalize(const CoordinateSequence& pts, std::size_t index)
{
    if (index >= pts.size() - 1 && pts.isRing()) {
        return 0;
    }
    return index;
}

CoordinateSequenceLocation
FacetSequence::nearestLocation(const CoordinateXY& p) const
{
    if (isPoint()) {
        // JTS uses index 0 / pts[0]; use start so a one-point facet
        // that is a subsequence of a longer CoordinateSequence is correct.
        return CoordinateSequenceLocation(pts, start, pts->getAt(start));
    }
    return nearestLocationOnLine(p);
}

CoordinateSequenceLocation
FacetSequence::nearestLocationOnLine(const CoordinateXY& pt) const
{
    if (clothoid) {
        const algorithm::clothoid::PlacedClothoid placed(*clothoid);
        const auto n = algorithm::clothoid::nearestToPoint(placed, pt);
        return CoordinateSequenceLocation(pts, start, Coordinate(n.onCurve));
    }
    if (arcs) {
        double minDistance = DoubleInfinity;
        std::size_t index = start;
        CoordinateXY nearestPt;
        for (std::size_t i = start; i + 2 < end; i += 2) {
            const CoordinateXY p = CircularArc(*pts, i).closestPoint(pt);
            const double dist = p.distance(pt);
            if (dist < minDistance) {
                minDistance = dist;
                nearestPt = p;
                index = i;
            }
        }
        return CoordinateSequenceLocation(pts, index, Coordinate(nearestPt));
    }

    double minDistance = DoubleInfinity;
    std::size_t index = start;
    Coordinate nearestPt;

    const Coordinate queryPt(pt);
    for (std::size_t i = start; i < end - 1; i++) {
        const Coordinate& q0 = pts->getAt(i);
        const Coordinate& q1 = pts->getAt(i + 1);
        double dist = Distance::pointToSegment(queryPt, q0, q1);
        if (dist < minDistance) {
            minDistance = dist;
            LineSegment seg(q0, q1);
            seg.closestPoint(queryPt, nearestPt);
            index = i;
            //-- segments are half-open, so 2nd endpoint belongs to next segment
            //-- except for last segment on non-closed sequence
            if (dist == 0.0 && queryPt.equals2D(q1)) {
                if (index < pts->size() - 1) {
                    index++;
                }
                //-- normalize index for a ring
                index = normalize(*pts, index);
            }
            if (minDistance <= 0.0) {
                break;
            }
        }
    }
    return CoordinateSequenceLocation(pts, index, nearestPt);
}

