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

#include <geos/algorithm/exactcurve/ExactCircularArc.h>
#include <geos/algorithm/exactcurve/ExactCurve.h>
#include <geos/constants.h>
#include <geos/geom/CircularString.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/LineString.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/util/IllegalArgumentException.h>

#include <cmath>
#include <optional>
#include <random>
#include <string>

using geos::algorithm::exactcurve::CookedArcs;
using geos::algorithm::exactcurve::ExactCircularArc;
using geos::algorithm::exactcurve::ExactCurve;
using geos::algorithm::exactcurve::Hit;
using geos::geom::CircularString;
using geos::geom::CoordinateXY;
using geos::geom::LineString;

namespace tut {

struct test_exact_circular_arc_data {
};

typedef test_group<test_exact_circular_arc_data> group;
typedef group::object object;

group test_exact_circular_arc_group(
    "geos::algorithm::exactcurve::ExactCircularArc");

// Semicircle (5,0)-(0,5)-(-5,0) length is 5π, not the chord 10.
template<>
template<>
void object::test<1>()
{
    set_test_name("semicircle length is five pi");
    ExactCircularArc a(CoordinateXY{5, 0}, CoordinateXY{0, 5}, CoordinateXY{-5, 0});
    ensure(a.isArc());
    ensure(std::fabs(a.radius() - 5.0) < 1.0e-12);
    ensure(std::fabs(a.sweep() - geos::MATH_PI) < 1.0e-12);
    ensure(std::fabs(a.length() - 5.0 * geos::MATH_PI) < 1.0e-12);
    ensure(a.chordLeArc());
    ensure(std::fabs(a.chordLength() - 10.0) < 1.0e-12);
}

template<>
template<>
void object::test<2>()
{
    set_test_name("collinear triple is a chord");
    ExactCircularArc a(CoordinateXY{0, 0}, CoordinateXY{1, 0}, CoordinateXY{3, 0});
    ensure(!a.isArc());
    ensure(std::fabs(a.length() - 3.0) < 1.0e-12);
    ensure(a.length() == a.chordLength());
    ensure(a.circularSegmentArea() == 0.0);
}

template<>
template<>
void object::test<3>()
{
    set_test_name("static length matches instance");
    CoordinateXY s{1, 0};
    CoordinateXY m{0, 1};
    CoordinateXY e{-1, 0};
    ensure(ExactCircularArc::lengthOf(s, m, e)
           == ExactCircularArc(s, m, e).length());
}

template<>
template<>
void object::test<4>()
{
    set_test_name("circular segment area of a half disc");
    ExactCircularArc a(CoordinateXY{5, 0}, CoordinateXY{0, 5}, CoordinateXY{-5, 0});
    ensure(std::fabs(a.circularSegmentArea() - 12.5 * geos::MATH_PI) < 1.0e-12);
}

template<>
template<>
void object::test<5>()
{
    set_test_name("arc-length centroid of a semicircle");
    ExactCircularArc a(CoordinateXY{5, 0}, CoordinateXY{0, 5}, CoordinateXY{-5, 0});
    const CoordinateXY c = a.arcLengthCentroid();
    ensure(std::fabs(c.x) < 1.0e-12);
    ensure(std::fabs(c.y - 10.0 / geos::MATH_PI) < 1.0e-12);
}

template<>
template<>
void object::test<6>()
{
    set_test_name("toLinear emits more than the chord");
    ExactCircularArc a(CoordinateXY{5, 0}, CoordinateXY{0, 5}, CoordinateXY{-5, 0});
    auto lin = a.toLinear(0.01);
    const auto* ls = dynamic_cast<const LineString*>(lin.get());
    ensure(ls != nullptr);
    ensure(ls->getNumPoints() > 2);
    ensure(ls->getCoordinateN(0).x == 5.0);
    ensure(ls->getCoordinateN(ls->getNumPoints() - 1).x == -5.0);
}

template<>
template<>
void object::test<7>()
{
    set_test_name("pointAt rejects out of range");
    ExactCircularArc a(CoordinateXY{5, 0}, CoordinateXY{0, 5}, CoordinateXY{-5, 0});
    try {
        a.pointAt(-0.1);
        fail("expected IAE");
    }
    catch (const geos::util::IllegalArgumentException&) {
    }
}

template<>
template<>
void object::test<8>()
{
    set_test_name("chord and inscribed polyline stay under exact length");
    std::mt19937 rng(2817130497u);
    std::uniform_real_distribution<double> box(-100.0, 100.0);
    int inscribedOver = 0;
    int chordOver = 0;
    constexpr int N = 8000;
    constexpr int nChord = 64;
    for (int i = 0; i < N; ++i) {
        ExactCircularArc a(CoordinateXY{box(rng), box(rng)},
                           CoordinateXY{box(rng), box(rng)},
                           CoordinateXY{box(rng), box(rng)});
        if (!a.chordLeArc()) {
            ++chordOver;
        }
        if (!a.isArc()) {
            continue;
        }
        const double exact = a.length();
        const double inscribed = nChord * 2.0 * a.radius()
            * std::sin(a.sweep() / (2.0 * nChord));
        if (inscribed > exact + 1.0e-12) {
            ++inscribedOver;
        }
    }
    ensure_equals(inscribedOver, 0);
    ensure_equals(chordOver, 0);
}

template<>
template<>
void object::test<9>()
{
    set_test_name("ExactCircularArc is-a ExactCurve");
    ExactCircularArc a(CoordinateXY{5, 0}, CoordinateXY{0, 5}, CoordinateXY{-5, 0});
    const ExactCurve& curve = a;
    const ExactCurve* p = &a;
    ensure(curve.isExact());
    ensure(p->isExact());
    ensure(curve.getStart().x == 5.0);
    ensure(curve.getEnd().x == -5.0);
    ensure(std::fabs(curve.length() - 5.0 * geos::MATH_PI) < 1.0e-12);
    const CoordinateXY mid = curve.pointAt(0.5);
    ensure(std::fabs(mid.x) < 1.0e-12);
    ensure(std::fabs(mid.y - 5.0) < 1.0e-12);
    auto lin = curve.toLinear(0.01);
    const auto* ls = dynamic_cast<const LineString*>(lin.get());
    ensure(ls != nullptr);
    ensure(ls->getNumPoints() > 2);

    ExactCircularArc chord(CoordinateXY{0, 0}, CoordinateXY{1, 0}, CoordinateXY{3, 0});
    const ExactCurve& chordCurve = chord;
    ensure(!chord.isArc());
    ensure(chordCurve.isExact());
    ensure(std::fabs(chordCurve.length() - 3.0) < 1.0e-12);
}

template<>
template<>
void object::test<10>()
{
    set_test_name("cook inserts one vertex on a crossing pair");
    geos::io::WKTReader reader;
    geos::io::WKTWriter writer;
    writer.setTrim(true);

    auto gA = reader.read<CircularString>("CIRCULARSTRING (-1 0, 0 1, 1 0)");
    auto gB = reader.read<CircularString>("CIRCULARSTRING (0 0, 1 1, 2 0)");
    ensure_equals(gA->getNumPoints(), static_cast<std::size_t>(3));
    ensure_equals(gB->getNumPoints(), static_cast<std::size_t>(3));
    ensure_equals(writer.write(*gA), std::string("CIRCULARSTRING (-1 0, 0 1, 1 0)"));
    ensure_equals(writer.write(*gB), std::string("CIRCULARSTRING (0 0, 1 1, 2 0)"));

    ExactCircularArc a(gA->getCoordinateN(0), gA->getCoordinateN(1), gA->getCoordinateN(2));
    ExactCircularArc b(gB->getCoordinateN(0), gB->getCoordinateN(1), gB->getCoordinateN(2));
    ensure(a.isArc());
    ensure(b.isArc());

    const auto h = ExactCircularArc::hit(a, b);
    ensure(h.has_value());
    ensure(std::fabs(h->point.x - 0.5) < 1.0e-10);
    ensure(std::fabs(h->point.y - std::sqrt(3.0) / 2.0) < 1.0e-10);

    const std::optional<CookedArcs> cooked = ExactCircularArc::cook(a, b);
    ensure(cooked.has_value());
    ensure(cooked->a0.isExact());
    ensure(cooked->a1.isExact());
    ensure(cooked->b0.isExact());
    ensure(cooked->b1.isExact());
    ensure(cooked->a0.getEnd().equals2D(cooked->hit.point));
    ensure(cooked->a1.getStart().equals2D(cooked->hit.point));
    ensure(cooked->b0.getEnd().equals2D(cooked->hit.point));
    ensure(cooked->b1.getStart().equals2D(cooked->hit.point));
    ensure(cooked->a0.isArc());
    ensure(cooked->a1.isArc());
    ensure(std::fabs(cooked->a0.center().x - a.center().x) < 1.0e-10);
    ensure(std::fabs(cooked->a0.center().y - a.center().y) < 1.0e-10);
    ensure(std::fabs(cooked->a1.center().x - a.center().x) < 1.0e-10);
    ensure(std::fabs(cooked->a1.center().y - a.center().y) < 1.0e-10);

    std::unique_ptr<CircularString> aOut;
    std::unique_ptr<CircularString> bOut;
    Hit inserted;
    ensure(ExactCircularArc::cook(*gA, *gB, aOut, bOut, inserted));
    ensure(aOut != nullptr);
    ensure(bOut != nullptr);
    ensure_equals("one hit vertex on A", aOut->getNumPoints(), static_cast<std::size_t>(5));
    ensure_equals("one hit vertex on B", bOut->getNumPoints(), static_cast<std::size_t>(5));
    ensure(aOut->hasCurvedComponents());
    ensure(bOut->hasCurvedComponents());
    ensure(std::fabs(inserted.point.x - 0.5) < 1.0e-10);
    ensure(std::fabs(inserted.point.y - std::sqrt(3.0) / 2.0) < 1.0e-10);
    ensure(aOut->getCoordinateN(2).equals2D(inserted.point));
    ensure(bOut->getCoordinateN(2).equals2D(inserted.point));
    ensure(aOut->getNumPoints() < 10);
    ensure(bOut->getNumPoints() < 10);
    ensure(writer.write(*aOut).find("CIRCULARSTRING") == 0);
    ensure(writer.write(*bOut).find("CIRCULARSTRING") == 0);

    ExactCircularArc a0(aOut->getCoordinateN(0), aOut->getCoordinateN(1), aOut->getCoordinateN(2));
    ExactCircularArc a1(aOut->getCoordinateN(2), aOut->getCoordinateN(3), aOut->getCoordinateN(4));
    ensure(a0.isExact());
    ensure(a1.isExact());
    ensure(a0.isArc());
    ensure(a1.isArc());
}

template<>
template<>
void object::test<11>()
{
    set_test_name("hit none when arcs miss");
    ExactCircularArc a(CoordinateXY{-1, 0}, CoordinateXY{0, 1}, CoordinateXY{1, 0});
    ExactCircularArc b(CoordinateXY{10, 0}, CoordinateXY{11, 1}, CoordinateXY{12, 0});
    ensure(!ExactCircularArc::hit(a, b).has_value());
    ensure(!ExactCircularArc::cook(a, b).has_value());
}

template<>
template<>
void object::test<12>()
{
    set_test_name("hit none when two proper intersections");
    // Nearly-full arcs of (0,0) r=1 and (1,0) r=1 meet at (1/2, ±√3/2).
    const double q = std::sqrt(2.0) / 2.0;
    ExactCircularArc a(CoordinateXY{1, 0}, CoordinateXY{-1, 0}, CoordinateXY{q, -q});
    ExactCircularArc b(CoordinateXY{2, 0}, CoordinateXY{0, 0}, CoordinateXY{1 + q, -q});
    ensure(a.isArc());
    ensure(b.isArc());
    ensure(!ExactCircularArc::hit(a, b).has_value());
    ensure(!ExactCircularArc::cook(a, b).has_value());
}

template<>
template<>
void object::test<13>()
{
    set_test_name("hit none on endpoint touch");
    ExactCircularArc a(CoordinateXY{-1, 0}, CoordinateXY{0, 1}, CoordinateXY{1, 0});
    ExactCircularArc b(CoordinateXY{1, 0}, CoordinateXY{2, 1}, CoordinateXY{3, 0});
    ensure(a.isArc());
    ensure(b.isArc());
    ensure(!ExactCircularArc::hit(a, b).has_value());
    ensure(!ExactCircularArc::cook(a, b).has_value());
}

template<>
template<>
void object::test<14>()
{
    set_test_name("hit none when cocircular");
    ExactCircularArc a(CoordinateXY{-1, 0}, CoordinateXY{0, 1}, CoordinateXY{1, 0});
    ExactCircularArc b(CoordinateXY{1, 0}, CoordinateXY{0, -1}, CoordinateXY{-1, 0});
    ensure(a.isArc());
    ensure(b.isArc());
    ensure(!ExactCircularArc::hit(a, b).has_value());
    ensure(!ExactCircularArc::cook(a, b).has_value());
}

template<>
template<>
void object::test<15>()
{
    set_test_name("hit none when a chord is collinear");
    ExactCircularArc chord(CoordinateXY{0, 0}, CoordinateXY{1, 0}, CoordinateXY{3, 0});
    ExactCircularArc arc(CoordinateXY{-1, 0}, CoordinateXY{0, 1}, CoordinateXY{1, 0});
    ensure(!chord.isArc());
    ensure(arc.isArc());
    ensure(!ExactCircularArc::hit(chord, arc).has_value());
    ensure(!ExactCircularArc::hit(arc, chord).has_value());
    ensure(!ExactCircularArc::cook(chord, arc).has_value());
}

} // namespace tut
