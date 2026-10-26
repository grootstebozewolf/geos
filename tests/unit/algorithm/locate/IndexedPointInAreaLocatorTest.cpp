#include <tut/tut.hpp>
// geos
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/algorithm/locate/IndexedPointInAreaLocator.h>
#include <geos/algorithm/locate/SimplePointInAreaLocator.h>
#include <geos/io/WKTReader.h>

// std
#include <memory>

using geos::geom::CoordinateXY;
using geos::geom::CoordinateSequence;
using geos::geom::GeometryFactory;
using geos::geom::Location;
using geos::algorithm::locate::IndexedPointInAreaLocator;
using geos::algorithm::locate::SimplePointInAreaLocator;

namespace tut {
//
// Test Group
//

// dummy data, not used
struct test_indexedpointinarealocator_data {
    const GeometryFactory& factory_ = *GeometryFactory::getDefaultInstance();
    geos::io::WKTReader reader_;

    void checkLocation(const std::string& wkt, const CoordinateXY& pt, Location expected)
    {
        auto g = reader_.read(wkt);
        IndexedPointInAreaLocator ipa(*g);
        ensure_equals(ipa.locate(&pt), expected);
        ensure_equals(SimplePointInAreaLocator::locate(pt, g.get()), expected);
    }
};

typedef test_group<test_indexedpointinarealocator_data> group;
typedef group::object object;

group test_indexedpointinarealocator_group("geos::algorithm::locate::IndexedPointInAreaLocator");

// Test all CoordinateSequence dimensions
template<>
template<>
void object::test<1>
()
{
    // XY
    CoordinateSequence seq_xy = CoordinateSequence::XY(0);
    seq_xy.add(0.0, 0.0);
    seq_xy.add(1.0, 0.0);
    seq_xy.add(1.0, 1.0);
    seq_xy.add(0.0, 1.0);
    seq_xy.add(0.0, 0.0);
    auto ls_xy = factory_.createLineString(seq_xy.clone());
    IndexedPointInAreaLocator ipa_xy(*ls_xy);

    CoordinateXY pt_boundary(0.5, 0);
    CoordinateXY pt_interior(0.5, 0.5);
    CoordinateXY pt_exterior(1.5, 0.5);

    ensure_equals(ipa_xy.locate(&pt_boundary), Location::BOUNDARY);
    ensure_equals(ipa_xy.locate(&pt_interior), Location::INTERIOR);
    ensure_equals(ipa_xy.locate(&pt_exterior), Location::EXTERIOR);

    // XYZ
    CoordinateSequence seq_xyz = CoordinateSequence::XYZ(0);
    seq_xyz.add(seq_xy);
    auto ls_xyz = factory_.createLineString(seq_xyz.clone());

    IndexedPointInAreaLocator ipa_xyz(*ls_xy);
    ensure_equals(ipa_xyz.locate(&pt_boundary), Location::BOUNDARY);
    ensure_equals(ipa_xyz.locate(&pt_interior), Location::INTERIOR);
    ensure_equals(ipa_xyz.locate(&pt_exterior), Location::EXTERIOR);

    // XYM
    CoordinateSequence seq_xym = CoordinateSequence::XYM(0);
    seq_xym.add(seq_xy);
    auto ls_xym = factory_.createLineString(seq_xym.clone());

    IndexedPointInAreaLocator ipa_xym(*ls_xy);
    ensure_equals(ipa_xym.locate(&pt_boundary), Location::BOUNDARY);
    ensure_equals(ipa_xym.locate(&pt_interior), Location::INTERIOR);
    ensure_equals(ipa_xym.locate(&pt_exterior), Location::EXTERIOR);

    // XYZM
    CoordinateSequence seq_xyzm = CoordinateSequence::XYZM(0);
    seq_xyzm.add(seq_xy);
    auto ls_xyzm = factory_.createLineString(seq_xyzm.clone());

    IndexedPointInAreaLocator ipa_xyzm(*ls_xy);
    ensure_equals(ipa_xyzm.locate(&pt_boundary), Location::BOUNDARY);
    ensure_equals(ipa_xyzm.locate(&pt_interior), Location::INTERIOR);
    ensure_equals(ipa_xyzm.locate(&pt_exterior), Location::EXTERIOR);
}

template<>
template<>
void object::test<2>
()
{
    set_test_name("CurvePolygon with a circular hole");

    const std::string wkt = "CURVEPOLYGON (COMPOUNDCURVE ((-10 100, 0 110, 10 100), "
                            "CIRCULARSTRING (10 100, 35 75, 10 50, 35 25, 10 0), (10 0, 0 -10, -10 0), "
                            "CIRCULARSTRING (-10 0, -35 25, -10 50, -35 75, -10 100)), "
                            "CIRCULARSTRING (-5 50, 0 55, 5 50, 0 45, -5 50))";

    checkLocation(wkt, CoordinateXY(30, 75), Location::INTERIOR);  // inside a bulging arc
    checkLocation(wkt, CoordinateXY(35, 75), Location::BOUNDARY);  // arc control point
    checkLocation(wkt, CoordinateXY(36, 75), Location::EXTERIOR);
    checkLocation(wkt, CoordinateXY(-12, 50), Location::EXTERIOR); // in the notch between two arcs
    checkLocation(wkt, CoordinateXY(0, 50), Location::EXTERIOR);   // inside the hole
    checkLocation(wkt, CoordinateXY(3, 54), Location::BOUNDARY);   // on the hole
}

template<>
template<>
void object::test<3>
()
{
    set_test_name("MultiSurface");

    const std::string wkt = "MULTISURFACE (CURVEPOLYGON (CIRCULARSTRING (-1 0, 0 1, 1 0, 0 -1, -1 0)), ((5 5, 6 5, 6 6, 5 5)))";

    checkLocation(wkt, CoordinateXY(0, 0.9), Location::INTERIOR);
    checkLocation(wkt, CoordinateXY(0.9, 0.9), Location::EXTERIOR);
    checkLocation(wkt, CoordinateXY(0, -1), Location::BOUNDARY);
    checkLocation(wkt, CoordinateXY(5.9, 5.5), Location::INTERIOR);
}

}

