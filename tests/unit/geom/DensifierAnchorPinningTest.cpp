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
#include <geos/io/WKTReader.h>

#include <cmath>
#include <memory>
#include <string>

using geos::algorithm::CurveToLineParams;
using geos::geom::CircularString;
using geos::geom::CoordinateXY;
using geos::geom::GeometryFactory;
using geos::io::WKTReader;

namespace tut {

struct test_densifier_anchor_data {
    GeometryFactory::Ptr factory_ = GeometryFactory::create();
    WKTReader wktreader_;

    std::unique_ptr<CircularString> readCS(const std::string& wkt)
    {
        return wktreader_.read<CircularString>(wkt);
    }

    bool containsExactly(const geos::geom::Geometry& g, const CoordinateXY& anchor)
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

// Witness: 15° steps on the unit semicircle put a computed vertex on
// the mid-control sweep. That vertex must be the supplied (0, 1), not
// cos(π/2) ≈ 6e-17.
template<>
template<>
void object::test<1>()
{
    set_test_name("mid control exact on sweep-angle tie");
    auto cs = readCS("CIRCULARSTRING (1 0, 0 1, -1 0)");
    auto ls = cs->getLinearized(CurveToLineParams::stepSizeDegrees(15));
    ensure("apex (0,1) survives getLinearized exactly",
           containsExactly(*ls, CoordinateXY{0, 1}));
}

template<>
template<>
void object::test<2>()
{
    set_test_name("mid control exact when it falls between vertices");
    auto cs = readCS("CIRCULARSTRING (1 0, 0 1, -1 0)");
    auto ls = cs->getLinearized(CurveToLineParams::stepSizeDegrees(16));
    ensure("apex (0,1) is inserted between samples",
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
