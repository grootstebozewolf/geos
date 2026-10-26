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

#pragma once

#include <geos/geom/Envelope.h>
#include <geos/geom/Coordinate.h>
#include <geos/operation/distance/CoordinateSequenceLocation.h>

#include <array>
#include <vector>

namespace geos {
namespace geom {
class Clothoid;
class CoordinateSequence;
}
}

namespace geos {
namespace operation {
namespace distance {

class FacetSequence {

private:

    const geom::CoordinateSequence* pts;
    const std::size_t start;
    const std::size_t end;
    // whether consecutive coordinates are circular arcs (p0, p1, p2) rather than segments
    const bool arcs;
    // a clothoid forming the single facet of this sequence, or null
    const geom::Clothoid* clothoid;

    /*
    * Unlike JTS, we store the envelope in the FacetSequence so
    * that it has a clear owner.  This is helpful when making a
    * tree of FacetSequence objects (FacetSequenceTreeBuilder)
    */
    geom::Envelope env;

    double computeDistanceLineLine(const FacetSequence& facetSeq,
                                   std::vector<geom::Coordinate> *locs) const;

    double computeDistancePointLine(const geom::Coordinate& pt,
                                    const FacetSequence& facetSeq,
                                    std::vector<geom::Coordinate> *locs) const;

    void updateNearestLocationsPointLine(const geom::Coordinate& pt,
                                         const geom::Coordinate& q0, const geom::Coordinate &q1,
                                         std::vector<geom::Coordinate> *locs) const;

    void updateNearestLocationsLineLine(const geom::Coordinate& p0, const geom::Coordinate& p1,
                                        const geom::Coordinate& q0, const geom::Coordinate &q1,
                                        std::vector<geom::Coordinate> *locs) const;

    CoordinateSequenceLocation nearestLocationOnLine(const geom::CoordinateXY& pt) const;

    std::size_t numElements() const;

    std::array<geom::CoordinateXY, 2> closestPoints(std::size_t i,
                                                    const FacetSequence& facetSeq, std::size_t j) const;

    geom::Envelope elementEnvelope(std::size_t i) const;

    double computeDistanceWithArcs(const FacetSequence& facetSeq,
                                   std::vector<geom::Coordinate>* locs) const;

    static std::size_t normalize(const geom::CoordinateSequence& pts, std::size_t index);

    void computeEnvelope();

public:

    const geom::Envelope* getEnvelope() const;

    const geom::Coordinate* getCoordinate(std::size_t index) const;

    std::size_t size() const;

    bool isPoint() const;

    double distance(const FacetSequence& facetSeq) const;

    /**
     * Creates a sequence of the facets between coordinates start and end (exclusive).
     *
     * @param pts the coordinates
     * @param start index of the first coordinate
     * @param end one past the index of the last coordinate
     * @param arcs whether the coordinates describe circular arcs, in which case
     *             start and end must lie on arc endpoints
     */
    FacetSequence(const geom::CoordinateSequence* pts, std::size_t start, std::size_t end, bool arcs = false);

    /**
     * Creates a sequence whose single facet is a clothoid.
     *
     * @param clothoid a non-empty clothoid
     */
    explicit FacetSequence(const geom::Clothoid* clothoid);

    /**
     * Computes the locations of the nearest points between this sequence
     * and another sequence.
     * The locations are presented in the same order as the input sequences.
     *
     * @return a pair of Coordinates for the nearest points
     */
    std::vector<geom::Coordinate> nearestLocations(const FacetSequence& facetSeq) const;

    /**
     * Computes the location of the nearest point of this sequence
     * to a query point.
     *
     * @param p the query point
     * @return the nearest location on this sequence
     */
    CoordinateSequenceLocation nearestLocation(const geom::CoordinateXY& p) const;

};

}
}
}
