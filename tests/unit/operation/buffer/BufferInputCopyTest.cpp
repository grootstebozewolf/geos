//
// Buffer must not mutate input coordinates.

#include <tut/tut.hpp>
#include <utility.h>

#include <geos/io/WKTReader.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>

using geos::geom::Coordinate;
using geos::geom::Geometry;
using geos::io::WKTReader;

namespace tut {
struct test_bufferinputcopy_data {
    WKTReader wktreader;
};
typedef test_group<test_bufferinputcopy_data> group;
typedef group::object object;

group test_bufferinputcopy_group("geos::operation::buffer::BufferInputCopy");

// testBufferDoesNotMutateInputCoordinates
template<>
template<>
void object::test<1>()
{
    auto input = wktreader.read("LINESTRING (0 0, 5 1, 10 0)");
    Coordinate orig1 = input->getCoordinates()->getAt(1);
    auto buf = input->buffer(1.0);
    (void) buf;
    Coordinate after = input->getCoordinates()->getAt(1);
    ensure(after.equals2D(orig1));
}

// testBufferThinLinearNoSpurious
template<>
template<>
void object::test<2>()
{
    auto thinRect = wktreader.read("POLYGON((0 0, 100 0, 100 0.1, 0 0.1, 0 0))");
    auto bufSmallNeg = thinRect->buffer(-0.01);
    ensure(bufSmallNeg->getNumGeometries() <= 1);
    if (!bufSmallNeg->isEmpty()) {
        ensure(bufSmallNeg->getArea() > 0);
    }
}

} // namespace tut
