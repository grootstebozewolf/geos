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
 * Last port: operation/distance/FacetSequenceTreeBuilder.java (f6187ee2 JTS-1.14)
 *
 **********************************************************************/

#include <geos/operation/distance/FacetSequenceTreeBuilder.h>
#include <geos/geom/Curve.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Point.h>
#include <geos/geom/SimpleCurve.h>

#include <algorithm>

using namespace geos::geom;
using namespace geos::index::strtree;

namespace geos {
namespace operation {
namespace distance {

std::unique_ptr<TemplateSTRtree<const FacetSequence*>>
FacetSequenceTreeBuilder::build(const Geometry* g)
{
    std::unique_ptr<TemplateSTRtree<const FacetSequence*>> tree(new FacetSequenceTree(computeFacetSequences(g)));

    tree->build();
    return tree;
}

std::vector<FacetSequence>
FacetSequenceTreeBuilder::computeFacetSequences(const Geometry* g)
{
    std::vector<FacetSequence> sections;

    class FacetSequenceAdder : public geom::GeometryComponentFilter {
        std::vector<FacetSequence>&  m_sections;

        public :
            FacetSequenceAdder(std::vector<FacetSequence> & p_sections) :
                m_sections(p_sections) {}
            void
            filter_ro(const Geometry* geom) override
            {
                if(const Curve* curve = dynamic_cast<const Curve*>(geom)) {
                    for (std::size_t i = 0; i < curve->getNumCurves(); i++) {
                        const SimpleCurve* section = curve->getCurveN(i);
                        const CoordinateSequence* seq = section->getCoordinatesRO();
                        if (section->getGeometryTypeId() == GEOS_CIRCULARSTRING) {
                            addArcFacetSequences(seq, m_sections);
                        }
                        else {
                            addFacetSequences(seq, m_sections);
                        }
                    }
                }
                else if(const Point* pt = dynamic_cast<const Point*>(geom)) {
                    const CoordinateSequence* seq = pt->getCoordinatesRO();
                    addFacetSequences(seq, m_sections);
                }
            }
    };

    FacetSequenceAdder facetSequenceAdder(sections);
    g->apply_ro(&facetSequenceAdder);

    return sections;
}

void
FacetSequenceTreeBuilder::addFacetSequences(
    const CoordinateSequence* pts,
    std::vector<FacetSequence>& sections)
{
    std::size_t i = 0;
    std::size_t size = pts->size();
    if ( size == 0 ) return;

    while(i <= size - 1) {
        std::size_t end = i + FACET_SEQUENCE_SIZE + 1;
        // if only one point remains after this section, include it in this
        // section
        if(end >= size - 1) {
            end = size;
        }
        sections.emplace_back(pts, i, end);
        i += FACET_SEQUENCE_SIZE;
    }
}

void
FacetSequenceTreeBuilder::addArcFacetSequences(
    const CoordinateSequence* pts,
    std::vector<FacetSequence>& sections)
{
    // Sections hold whole arcs, so they start and end on arc endpoints.
    const std::size_t size = pts->size();
    if (size < 3) return;

    const std::size_t step = 2 * (FACET_SEQUENCE_SIZE / 2);
    for (std::size_t i = 0; i + 1 < size; i += step) {
        const std::size_t end = std::min(i + step + 1, size);
        sections.emplace_back(pts, i, end, true);
    }
}

}
}
}

