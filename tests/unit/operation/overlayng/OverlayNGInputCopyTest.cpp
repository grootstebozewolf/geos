//
// OverlayNG on linear inputs must not mutate input coordinates.

#include <tut/tut.hpp>
#include <utility.h>

#include <geos/io/WKTReader.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Geometry.h>
#include <geos/operation/overlayng/OverlayNG.h>

using geos::geom::Coordinate;
using geos::geom::Geometry;
using geos::io::WKTReader;
using geos::operation::overlayng::OverlayNG;

namespace tut {
struct test_overlaynginputcopy_data {
    WKTReader wktreader;
};
typedef test_group<test_overlaynginputcopy_data> group;
typedef group::object object;

group test_overlaynginputcopy_group("geos::operation::overlayng::OverlayNGInputCopy");

// testOverlayNGDoesNotMutateInputCoordinatesLinear
template<>
template<>
void object::test<1>()
{
    auto inputA = wktreader.read("LINESTRING (0 0, 10 0, 20 0)");
    auto inputB = wktreader.read("LINESTRING (5 -5, 5 5)");
    Coordinate orig = inputA->getCoordinates()->getAt(1);
    auto result = OverlayNG::overlay(inputA.get(), inputB.get(), OverlayNG::INTERSECTION);
    (void) result;
    Coordinate after = inputA->getCoordinates()->getAt(1);
    ensure(after.equals2D(orig));
}

} // namespace tut
