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

#include <geos/algorithm/distance/DiscreteHausdorffDistance.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/io/WKTReader.h>

#include <cmath>
#include <memory>
#include <string>

using geos::algorithm::distance::DiscreteHausdorffDistance;
using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::geom::PrecisionModel;
using geos::io::WKTReader;

namespace tut {

struct test_dhd_curve_data {
    PrecisionModel pm;
    GeometryFactory::Ptr gf;
    WKTReader reader;

    test_dhd_curve_data()
        : pm(), gf(GeometryFactory::create(&pm)), reader(gf.get())
    {}

    double oriented(const std::string& wktA, const std::string& wktB)
    {
        std::unique_ptr<Geometry> a(reader.read(wktA));
        std::unique_ptr<Geometry> b(reader.read(wktB));
        DiscreteHausdorffDistance dhd(*a, *b);
        return dhd.orientedDistance();
    }
};

typedef test_group<test_dhd_curve_data> group;
typedef group::object object;

group test_dhd_curve_group(
    "geos::algorithm::distance::DiscreteHausdorffDistance::curve");

// Witness: CIRCULARSTRING (0 0, 2 3, 10 0) → LINESTRING (0 0, 10 0)
// is the apex √949/6 − 7/6, not the far-end chord 10 or mid-control 3.
template<>
template<>
void object::test<1>()
{
    set_test_name("arc to segment is apex not far end");
    const double apex = std::sqrt(949.0) / 6.0 - 7.0 / 6.0;
    const double got = oriented(
        "CIRCULARSTRING (0 0, 2 3, 10 0)",
        "LINESTRING (0 0, 10 0)");
    ensure(std::fabs(got - apex) < 1.0e-9);
    ensure(got < 9.0);
}

template<>
template<>
void object::test<2>()
{
    set_test_name("two discs match circle-to-circle");
    const double got = oriented(
        "CURVEPOLYGON (CIRCULARSTRING (5 0, 0 5, -5 0, 0 -5, 5 0))",
        "CURVEPOLYGON (CIRCULARSTRING (12 0, 7 5, 2 0, 7 -5, 12 0))");
    // |d + r1 − r2| for r=5 discs 7 apart: directed HD is 7, not 10.
    ensure(std::fabs(got - 7.0) < 1.0e-9);
}

} // namespace tut
