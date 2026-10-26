//
// Test Suite for geos::algorithm::PointLocator
// Ported from JTS junit/algorithm/PointLocator.java

#include <tut/tut.hpp>
// geos
#include <geos/io/WKTReader.h>
#include <geos/algorithm/locate/SimplePointInAreaLocator.h>
#include <geos/algorithm/PointLocator.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h> // required for use in unique_ptr
#include <geos/geom/Coordinate.h>
// std
#include <string>
#include <memory>

namespace geos {
namespace geom {
class Geometry;
}
}

using namespace geos::geom; // for Location

namespace tut {
//
// Test Group
//

// dummy data, not used
struct test_pointlocator_data {};

typedef test_group<test_pointlocator_data> group;
typedef group::object object;

group test_pointlocator_group("geos::algorithm::PointLocator");

// These are static to avoid namespace pollution
// The struct test_*_data above is probably there
// for the same reason...
//
static PrecisionModel pm;
static GeometryFactory::Ptr gf = GeometryFactory::create(&pm);
static geos::io::WKTReader reader(gf.get());

typedef std::unique_ptr<Geometry> GeomPtr;

void
runPtLocator(Location expected, const Coordinate& pt,
             const std::string& wkt)
{
    GeomPtr geom(reader.read(wkt));
    geos::algorithm::PointLocator pointLocator;
    Location loc = pointLocator.locate(pt, geom.get());
    ensure_equals(loc, expected);
}

void
runSimplePtLocator(Location expected, const Coordinate& pt,
             const std::string& wkt)
{
    GeomPtr geom(reader.read(wkt));
    Location loc = geos::algorithm::locate::SimplePointInAreaLocator::locate(pt, geom.get());
    ensure_equals(loc, expected);
}

//
// Test Cases
//

// 1 - Test box
template<>
template<>
void object::test<1>
()
{
    runPtLocator(Location::INTERIOR, Coordinate(10, 10),
                 "POLYGON ((0 0, 0 20, 20 20, 20 0, 0 0))");
}

// 2 - Test complex ring
template<>
template<>
void object::test<2>
()
{
    runPtLocator(Location::INTERIOR, Coordinate(0, 0),
                 "POLYGON ((-40 80, -40 -80, 20 0, 20 -100, 40 40, 80 -80, 100 80, 140 -20, 120 140, 40 180,     60 40, 0 120, -20 -20, -40 80))");
}

// 3 - Test PointLocator LinearRing LineString
template<>
template<>
void object::test<3>
()
{
    runPtLocator(Location::BOUNDARY, Coordinate(0, 0),
                 "GEOMETRYCOLLECTION( LINESTRING(0 0, 10 10), LINEARRING(10 10, 10 20, 20 10, 10 10))");
}

// 4 - Test PointLocator Point inside LinearRing
template<>
template<>
void object::test<4>
()
{
    runPtLocator(Location::EXTERIOR, Coordinate(11, 11),
                 "LINEARRING(10 10, 10 20, 20 10, 10 10)");
}

// 5 - TestPointLocator Point inside MultiPoint
template<>
template<>
void object::test<5>
()
{
    runPtLocator(Location::INTERIOR, Coordinate(0, 0),
                 "MULTIPOINT ((1 1), (0 0))");
}

// 6 - TestPointLocator Point inside GeometryCollection containing MultiPolygon, using SimplePointInAreaLocator.
template<>
template<>
void object::test<6>
()
{
    runSimplePtLocator(Location::INTERIOR, Coordinate(0, 0),
                       "GEOMETRYCOLLECTION (MULTIPOLYGON (((-1 -1, 1 -1, 1 1, -1 1, -1 -1))))");
}

// 7 - Endpoints of an open CircularString are its boundary
template<>
template<>
void object::test<7>
()
{
    runPtLocator(Location::BOUNDARY, Coordinate(-1, 0), "CIRCULARSTRING (-1 0, 0 1, 1 0)");
    runPtLocator(Location::BOUNDARY, Coordinate(1, 0), "CIRCULARSTRING (-1 0, 0 1, 1 0)");
    runPtLocator(Location::INTERIOR, Coordinate(0, 1), "CIRCULARSTRING (-1 0, 0 1, 1 0)");
    runPtLocator(Location::EXTERIOR, Coordinate(0, 0), "CIRCULARSTRING (-1 0, 0 1, 1 0)");
}

// 8 - A closed CircularString has no boundary
template<>
template<>
void object::test<8>
()
{
    runPtLocator(Location::INTERIOR, Coordinate(-1, 0), "CIRCULARSTRING (-1 0, 0 1, 1 0, 0 -1, -1 0)");
    runPtLocator(Location::EXTERIOR, Coordinate(0, 0), "CIRCULARSTRING (-1 0, 0 1, 1 0, 0 -1, -1 0)");
}

// 9 - Endpoints of an open CompoundCurve are its boundary; a closed one has none
template<>
template<>
void object::test<9>
()
{
    runPtLocator(Location::BOUNDARY, Coordinate(-1, 0), "COMPOUNDCURVE (CIRCULARSTRING (-1 0, 0 1, 1 0), (1 0, 2 0))");
    runPtLocator(Location::BOUNDARY, Coordinate(2, 0), "COMPOUNDCURVE (CIRCULARSTRING (-1 0, 0 1, 1 0), (1 0, 2 0))");
    runPtLocator(Location::INTERIOR, Coordinate(1, 0), "COMPOUNDCURVE (CIRCULARSTRING (-1 0, 0 1, 1 0), (1 0, 2 0))");
    runPtLocator(Location::INTERIOR, Coordinate(-1, 0), "COMPOUNDCURVE (CIRCULARSTRING (-1 0, 0 1, 1 0), (1 0, -1 0))");
}

// 10 - MultiCurve follows the mod-2 boundary rule
template<>
template<>
void object::test<10>
()
{
    const std::string wkt = "MULTICURVE (CIRCULARSTRING (-1 0, 0 1, 1 0), (1 0, 2 0))";
    runPtLocator(Location::INTERIOR, Coordinate(1, 0), wkt);
    runPtLocator(Location::BOUNDARY, Coordinate(2, 0), wkt);
    runPtLocator(Location::INTERIOR, Coordinate(0, 1), wkt);
    runPtLocator(Location::EXTERIOR, Coordinate(0, 0), wkt);
}

// 11 - MultiSurface and CurvePolygon elements of a collection
template<>
template<>
void object::test<11>
()
{
    const std::string ms = "MULTISURFACE (CURVEPOLYGON (CIRCULARSTRING (-1 0, 0 1, 1 0, 0 -1, -1 0)), "
                           "((5 5, 6 5, 6 6, 5 5)))";
    runPtLocator(Location::INTERIOR, Coordinate(0, 0.9), ms);
    runPtLocator(Location::BOUNDARY, Coordinate(0, 1), ms);
    runPtLocator(Location::EXTERIOR, Coordinate(0.9, 0.9), ms);

    runPtLocator(Location::INTERIOR, Coordinate(0, 0.5),
                 "GEOMETRYCOLLECTION (CURVEPOLYGON (COMPOUNDCURVE (CIRCULARSTRING (-1 0, 0 1, 1 0), (1 0, -1 0))), "
                 "POINT (5 5))");
}

} // namespace tut

