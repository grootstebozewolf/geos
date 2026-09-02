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

#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>

#include <geos/algorithm/CurveToLineParams.h>
#include <geos/geom/CircularString.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/io/WKTReader.h>

#include <memory>
#include <string>

using geos::algorithm::CurveToLineParams;
using geos::geom::CircularString;
using geos::geom::CoordinateXY;
using geos::geom::GeometryFactory;
using geos::geom::LineString;
using geos::io::WKTReader;

namespace tut {

// topic: arc
// Port of JTS f6347444.
struct test_densifier_anchor_data {
    GeometryFactory::Ptr factory_ = GeometryFactory::create();
    WKTReader wktreader_;

    std::unique_ptr<CircularString> readCS(const std::string& wkt)
    {
        return wktreader_.read<CircularString>(wkt);
    }

    bool containsExactly(const geos::geom::LineString& g, const CoordinateXY& anchor)
    {
        const auto* seq = g.getCoordinatesRO();
        for (std::size_t i = 0; i < seq->size(); ++i) {
            const CoordinateXY& p = seq->getAt<CoordinateXY>(i);
            if (p.equals2D(anchor)) {
                return true;
            }
        }
        return false;
    }
};

typedef test_group<test_densifier_anchor_data> group;
typedef group::object object;

group test_densifier_anchor_group(
    "geos::geom::CircularString::getLinearized::anchors");

// Witness: CIRCULARSTRING (1 0, 0 1, -1 0) after getLinearized(0.01)
// keeps exact (0, 1), not cos(π/2).
template<>
template<>
void object::test<1>()
{
    set_test_name("mid control exact at tolerance 0.01");
    auto cs = readCS("CIRCULARSTRING (1 0, 0 1, -1 0)");
    auto ls = cs->getLinearized(CurveToLineParams::maxDeviation(0.01));
    ensure("apex (0,1) survives getLinearized exactly",
           containsExactly(*ls, CoordinateXY{0, 1}));
}

template<>
template<>
void object::test<2>()
{
    set_test_name("mid control exact at tolerance 0.012");
    auto cs = readCS("CIRCULARSTRING (1 0, 0 1, -1 0)");
    auto ls = cs->getLinearized(CurveToLineParams::maxDeviation(0.012));
    ensure("apex (0,1) is present between samples",
           containsExactly(*ls, CoordinateXY{0, 1}));
}

template<>
template<>
void object::test<3>()
{
    set_test_name("all controls of a two-arc string survive");
    auto cs = readCS("CIRCULARSTRING (1 0, 0 1, -1 0, -2 -1, -3 0)");
    auto ls = cs->getLinearized(CurveToLineParams::maxDeviation(0.01));
    ensure(containsExactly(*ls, CoordinateXY{1, 0}));
    ensure(containsExactly(*ls, CoordinateXY{0, 1}));
    ensure(containsExactly(*ls, CoordinateXY{-1, 0}));
    ensure(containsExactly(*ls, CoordinateXY{-2, -1}));
    ensure(containsExactly(*ls, CoordinateXY{-3, 0}));
}

} // namespace tut
