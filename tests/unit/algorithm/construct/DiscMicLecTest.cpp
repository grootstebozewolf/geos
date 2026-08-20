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

#include <geos/algorithm/construct/LargestEmptyCircle.h>
#include <geos/algorithm/construct/MaximumInscribedCircle.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/io/WKTReader.h>

#include <cmath>
#include <memory>
#include <string>

using geos::algorithm::construct::LargestEmptyCircle;
using geos::algorithm::construct::MaximumInscribedCircle;
using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::geom::PrecisionModel;
using geos::io::WKTReader;

namespace tut {

// topic: arc
// Port of JTS f24cb33d.
struct test_disc_mic_lec_data {
    PrecisionModel pm;
    GeometryFactory::Ptr gf;
    WKTReader reader;

    test_disc_mic_lec_data()
        : pm(), gf(GeometryFactory::create(&pm)), reader(gf.get())
    {}
};

typedef test_group<test_disc_mic_lec_data> group;
typedef group::object object;

group test_disc_mic_lec_group(
    "geos::algorithm::construct::MaximumInscribedCircle::disc");

// Witness: CURVEPOLYGON (CIRCULARSTRING (5 0, 0 5, -5 0, 0 -5, 5 0))
// MIC/LEC radius is 5.0, not 5/√2.
template<>
template<>
void object::test<1>()
{
    set_test_name("radius-5 disc MIC is 5 not 5/sqrt(2)");
    std::unique_ptr<Geometry> disc(reader.read(
        "CURVEPOLYGON (CIRCULARSTRING (5 0, 0 5, -5 0, 0 -5, 5 0))"));
    auto radius = MaximumInscribedCircle::getRadiusLine(disc.get(), 0.01);
    ensure(std::fabs(radius->getLength() - 5.0) < 1.0e-9);
    ensure(radius->getLength() > 4.0);
}

template<>
template<>
void object::test<2>()
{
    set_test_name("radius-5 disc LEC is 5");
    std::unique_ptr<Geometry> disc(reader.read(
        "CURVEPOLYGON (CIRCULARSTRING (5 0, 0 5, -5 0, 0 -5, 5 0))"));
    auto radius = LargestEmptyCircle::getRadiusLine(disc.get(), 0.01);
    ensure(std::fabs(radius->getLength() - 5.0) < 1.0e-9);
}

} // namespace tut
