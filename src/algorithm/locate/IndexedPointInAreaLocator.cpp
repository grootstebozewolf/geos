/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 * Copyright (C) 2018 Daniel Baston <dbaston@gmail.com>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/


#include <geos/algorithm/locate/IndexedPointInAreaLocator.h>
#include <geos/algorithm/CircularArcs.h>
#include <geos/geom/CircularArc.h>
#include <geos/geom/Clothoid.h>
#include <geos/geom/Curve.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/SimpleCurve.h>
#include <geos/geom/Surface.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/LineString.h>
#include <geos/geom/LineSegment.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/index/strtree/Interval.h>
#include <geos/index/strtree/TemplateSTRtree.h>
#include <geos/util.h>
#include <geos/algorithm/RayCrossingCounter.h>
#include <geos/index/ItemVisitor.h>
#include <algorithm/ClothoidGeometry.h>

#include <algorithm>

using geos::geom::CoordinateXY;

namespace geos {
namespace algorithm {
namespace locate {
//
// private:
//
IndexedPointInAreaLocator::IntervalIndexedGeometry::IntervalIndexedGeometry(const geom::Geometry& g)
{
    init(g);
}

namespace {

// Collect the closed curves of a geometry: the rings of its surfaces and any
// closed curve (such as a LinearRing) appearing on its own.
void
getClosedCurves(const geom::Geometry& g, std::vector<const geom::Curve*>& rings)
{
    if (g.isEmpty()) {
        return;
    }

    if (g.isCollection()) {
        for (std::size_t i = 0; i < g.getNumGeometries(); i++) {
            getClosedCurves(*g.getGeometryN(i), rings);
        }
        return;
    }

    if (g.getDimension() == geom::Dimension::A) {
        const geom::Surface& surface = static_cast<const geom::Surface&>(g);
        getClosedCurves(*surface.getExteriorRing(), rings);
        for (std::size_t i = 0; i < surface.getNumInteriorRing(); i++) {
            getClosedCurves(*surface.getInteriorRingN(i), rings);
        }
        return;
    }

    if (g.getDimension() == geom::Dimension::L) {
        const geom::Curve& curve = static_cast<const geom::Curve&>(g);
        //-- only include rings of Polygons or LinearRings
        if (curve.isClosed()) {
            rings.push_back(&curve);
        }
    }
}

}

void
IndexedPointInAreaLocator::IntervalIndexedGeometry::init(const geom::Geometry& g)
{
    std::vector<const geom::Curve*> rings;
    getClosedCurves(g, rings);

    // pre-compute size of segment, arc and clothoid vectors
    std::size_t nsegs = 0;
    std::size_t narcs = 0;
    std::size_t nclothoids = 0;
    for (const geom::Curve* ring : rings) {
        for (std::size_t i = 0; i < ring->getNumCurves(); i++) {
            const geom::SimpleCurve* section = ring->getCurveN(i);
            const std::size_t npts = section->getCoordinatesRO()->size();
            if (section->getGeometryTypeId() == geom::GEOS_CIRCULARSTRING) {
                narcs += (npts - 1) / 2;
            }
            else if (section->getGeometryTypeId() == geom::GEOS_CLOTHOID) {
                nclothoids++;
            }
            else {
                nsegs += npts - 1;
            }
        }
    }
    index = decltype(index)(10, nsegs);
    arcIndex = decltype(arcIndex)(10, narcs);
    clothoidIndex = decltype(clothoidIndex)(10, nclothoids);

    for (const geom::Curve* ring : rings) {
        for (std::size_t i = 0; i < ring->getNumCurves(); i++) {
            const geom::SimpleCurve* section = ring->getCurveN(i);
            if (section->getGeometryTypeId() == geom::GEOS_CIRCULARSTRING) {
                addArcs(section->getCoordinatesRO());
            }
            else if (section->getGeometryTypeId() == geom::GEOS_CLOTHOID) {
                addClothoid(static_cast<const geom::Clothoid*>(section));
            }
            else {
                addLine(section->getCoordinatesRO());
            }
        }
    }
}

void
IndexedPointInAreaLocator::IntervalIndexedGeometry::addLine(const geom::CoordinateSequence* pts)
{
    for(std::size_t i = 1, ni = pts->size(); i < ni; i++) {
        SegmentView seg(&pts->getAt<CoordinateXY>(i-1), &pts->getAt<CoordinateXY>(i));
        auto r = std::minmax(seg.p0().y, seg.p1().y);

        index.insert(index::strtree::Interval(r.first, r.second), seg);
    }
}

void
IndexedPointInAreaLocator::IntervalIndexedGeometry::addArcs(const geom::CoordinateSequence* pts)
{
    for (std::size_t i = 2, ni = pts->size(); i < ni; i += 2) {
        geom::Envelope env;
        CircularArcs::expandEnvelope(env,
                                     pts->getAt<CoordinateXY>(i-2),
                                     pts->getAt<CoordinateXY>(i-1),
                                     pts->getAt<CoordinateXY>(i));

        arcIndex.insert(index::strtree::Interval(env.getMinY(), env.getMaxY()), ArcView{pts, i-2});
    }
}

void
IndexedPointInAreaLocator::IntervalIndexedGeometry::addClothoid(const geom::Clothoid* clothoid)
{
    if (clothoid->isEmpty()) {
        return;
    }
    // A point within the evaluation tolerance of the clothoid is on it.
    const double tol = algorithm::clothoid::PlacedClothoid(*clothoid).tolerance();
    const geom::Envelope* env = clothoid->getEnvelopeInternal();
    clothoidIndex.insert(index::strtree::Interval(env->getMinY() - tol, env->getMaxY() + tol), clothoid);
}


void
IndexedPointInAreaLocator::buildIndex(const geom::Geometry& g)
{
    index = detail::make_unique<IntervalIndexedGeometry>(g);
}


//
// protected:
//

//
// public:
//
IndexedPointInAreaLocator::IndexedPointInAreaLocator(const geom::Geometry& g)
    :	areaGeom(g)
{
}

geom::Location
IndexedPointInAreaLocator::locate(const geom::CoordinateXY* /*const*/ p)
{
    if (index == nullptr) {
        buildIndex(areaGeom);
    }

    algorithm::RayCrossingCounter rcc(*p);

    index->query(p->y, p->y, [&rcc](const SegmentView& ls) {
        rcc.countSegment(ls.p0(), ls.p1());
    });

    index->queryArcs(p->y, p->y, [&rcc](const ArcView& av) {
        rcc.countArc(geom::CircularArc(*av.seq, av.pos));
    });

    index->queryClothoids(p->y, p->y, [&rcc](const geom::Clothoid* clothoid) {
        rcc.countClothoid(*clothoid);
    });

    return rcc.getLocation();
}


} // geos::algorithm::locate
} // geos::algorithm
} // geos
