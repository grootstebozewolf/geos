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

#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/IntersectionMatrix.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/io/WKTReader.h>

#include <memory>
#include <string>

using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::geom::PrecisionModel;
using geos::io::WKTReader;

namespace tut {

struct test_disc_touch_data {
    PrecisionModel pm;
    GeometryFactory::Ptr gf;
    WKTReader reader;

    test_disc_touch_data()
        : pm(), gf(GeometryFactory::create(&pm)), reader(gf.get())
    {}
};

typedef test_group<test_disc_touch_data> group;
typedef group::object object;

group test_disc_touch_group("geos::geom::Geometry::touches::disc");

// Witness: r=5 discs with 3-4-5 centres kiss at (4, 3); DE-9IM FF2F01212.
template<>
template<>
void object::test<1>()
{
    set_test_name("external tangent discs are FF2F01212");
    std::unique_ptr<Geometry> a(reader.read(
        "CURVEPOLYGON (CIRCULARSTRING (-5 0, 0 5, 5 0, 0 -5, -5 0))"));
    std::unique_ptr<Geometry> b(reader.read(
        "CURVEPOLYGON (CIRCULARSTRING (13 6, 8 11, 3 6, 8 1, 13 6))"));
    auto im = a->relate(b.get());
    ensure_equals(im->toString(), std::string("FF2F01212"));
    ensure(a->touches(b.get()));
    ensure(b->touches(a.get()));
}

} // namespace tut
