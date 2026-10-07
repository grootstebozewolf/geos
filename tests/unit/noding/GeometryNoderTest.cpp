#include <tut/tut.hpp>
#include <cmath>

#include <geos/geom/Geometry.h>
#include <geos/geom/SimpleCurve.h>
#include <geos/io/WKTWriter.h>
#include <geos/noding/GeometryNoder.h>

#include "utility.h"

using geos::geom::Geometry;
using geos::noding::GeometryNoder;

namespace tut {

struct test_geometrynoder_data {
    geos::io::WKTReader reader_;

    // Node a and b, then check the piece types and that each hit
    // point is an endpoint of `ends` pieces.
    void checkArcNoding(const std::string& a, const std::string& b, bool preserve,
                        std::size_t nArcs, std::size_t nLines,
                        const std::vector<geos::geom::CoordinateXY>& hits, std::size_t ends)
    {
        auto g1 = reader_.read(a);
        auto g2 = reader_.read(b);
        GeometryNoder noder(*g1, *g2);
        noder.setPreserveCompoundCurves(preserve);
        auto result = noder.getNoded();
        ensure(result != nullptr);

        std::size_t arcs = 0, lines = 0;
        for (std::size_t i = 0; i < result->getNumGeometries(); i++) {
            auto t = result->getGeometryN(i)->getGeometryTypeId();
            if (t == geos::geom::GEOS_CIRCULARSTRING) arcs++;
            else if (t == geos::geom::GEOS_LINESTRING) lines++;
        }
        ensure_equals("circular pieces: " + result->toString(), arcs, nArcs);
        ensure_equals("linear pieces: " + result->toString(), lines, nLines);

        for (const auto& hit : hits) {
            std::size_t found = 0;
            for (std::size_t i = 0; i < result->getNumGeometries(); i++) {
                auto c = static_cast<const geos::geom::SimpleCurve*>(result->getGeometryN(i));
                if (c->getCoordinatesRO()->front<geos::geom::CoordinateXY>().distance(hit) < 1e-9) found++;
                if (c->getCoordinatesRO()->back<geos::geom::CoordinateXY>().distance(hit) < 1e-9) found++;
            }
            ensure_equals("pieces ending at hit: " + result->toString(), found, ends);
        }
    }
};

typedef test_group<test_geometrynoder_data> group;
typedef group::object object;

group test_geometrynoder_group("geos::noding::GeometryNoder");

template<>
template<>
void object::test<1>()
{
    set_test_name("single input") ;

    auto input = reader_.read("MULTILINESTRING ((0 0, 10 10, 10 0, 0 10))");

    auto result = GeometryNoder::node(*input);
    ensure(result != nullptr);

    auto expected = reader_.read("MULTILINESTRING ((0 0, 5 5), (5 5, 10 10, 10 0, 5 5), (5 5, 0 10))");

    ensure_equals_exact_geometry_xyzm( result.get(), expected.get(), 0);
}

template<>
template<>
void object::test<2>()
{
    set_test_name("two inputs, output all edges") ;

    auto input1 = reader_.read("LINESTRING (0 0, 10 10)");
    auto input2 = reader_.read("LINESTRING (0 10, 10 0)");

    auto result = GeometryNoder::node(*input1, *input2);
    ensure(result != nullptr);

    auto expected = reader_.read("MULTILINESTRING ((0 0, 5 5), (5 5, 10 10), (0 10, 5 5), (5 5, 10 0))");

    ensure_equals_exact_geometry_xyzm(result.get(), expected.get(), 0);
}

template<>
template<>
void object::test<3>()
{
    set_test_name("two inputs, only output edges from first") ;

    auto input1 = reader_.read("LINESTRING (0 0, 10 10)");
    auto input2 = reader_.read("LINESTRING (0 10, 10 0)");

    GeometryNoder noder(*input1, *input2);
    noder.setOnlyFirstGeomEdges(true);

    auto result = noder.getNoded();
    ensure(result != nullptr);

    auto expected = reader_.read("MULTILINESTRING ((0 0, 5 5), (5 5, 10 10))");

    ensure_equals_exact_geometry_xyzm(result.get(), expected.get(), 0);
}

// Upper half circles centred at (0 0) and (5 0) cross once at (2.5, 5*sqrt(3)/2).
template<>
template<>
void object::test<4>()
{
    set_test_name("two crossing circular strings stay circular");

    for (bool preserve : {false, true}) {
        checkArcNoding("CIRCULARSTRING (-5 0, 0 5, 5 0)",
                       "CIRCULARSTRING (0 0, 5 5, 10 0)",
                       preserve, 4, 0, {{2.5, 2.5 * std::sqrt(3.0)}}, 4);
    }
}

// The line x = 3 hits the arc at (3 4).
template<>
template<>
void object::test<5>()
{
    set_test_name("circular string crossing a line string");

    for (bool preserve : {false, true}) {
        checkArcNoding("CIRCULARSTRING (-5 0, 0 5, 5 0)",
                       "LINESTRING (3 -1, 3 10)",
                       preserve, 2, 2, {{3, 4}}, 4);
    }
}

// A full circle cut by a vertical line through its centre.
// The start point (-3 0) stays a node, so the circle gives three arcs.
template<>
template<>
void object::test<6>()
{
    set_test_name("full circle noded by a line");

    for (bool preserve : {false, true}) {
        checkArcNoding("CIRCULARSTRING (-3 0, 3 0, -3 0)",
                       "LINESTRING (0 -5, 0 5)",
                       preserve, 3, 3, {{0, 3}, {0, -3}}, 4);
    }
}

} // namespace tut
