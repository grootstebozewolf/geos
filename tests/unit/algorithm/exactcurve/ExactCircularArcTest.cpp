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
#include <geos/geom/Coordinate.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/LineString.h>
#include <geos/util/IllegalArgumentException.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using geos::algorithm::exactcurve::ExactCircularArc;
using geos::algorithm::exactcurve::ExactCurve;
using geos::geom::CoordinateXY;
using geos::geom::LineString;

namespace tut {

// topic: arc
// Port of JTS 9797c2c4.
struct test_exact_circular_arc_data {
};

typedef test_group<test_exact_circular_arc_data> group;
typedef group::object object;

group test_exact_circular_arc_group(
    "geos::algorithm::exactcurve::ExactCircularArc");

// Witness: semicircle (5,0)-(0,5)-(-5,0) length is 5π, not the chord 10.
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
    set_test_name("colinear is chord");
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
    set_test_name("circular segment area half disc");
    ExactCircularArc a(CoordinateXY{5, 0}, CoordinateXY{0, 5}, CoordinateXY{-5, 0});
    ensure(std::fabs(a.circularSegmentArea() - 12.5 * geos::MATH_PI) < 1.0e-12);
}

template<>
template<>
void object::test<5>()
{
    set_test_name("arc-length centroid semicircle");
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
    set_test_name("L1 L2 hard zero on random windows");
    std::mt19937 rng(2817130497u);
    std::uniform_real_distribution<double> box(-100.0, 100.0);
    int l1Hard = 0;
    int l2Hard = 0;
    constexpr int N = 8000;
    constexpr int nChord = 64;
    for (int i = 0; i < N; ++i) {
        ExactCircularArc a(CoordinateXY{box(rng), box(rng)},
                           CoordinateXY{box(rng), box(rng)},
                           CoordinateXY{box(rng), box(rng)});
        if (!a.chordLeArc()) {
            ++l2Hard;
        }
        if (!a.isArc()) {
            continue;
        }
        const double exact = a.length();
        const double inscribed = nChord * 2.0 * a.radius()
            * std::sin(a.sweep() / (2.0 * nChord));
        if (inscribed > exact + 1.0e-12) {
            ++l1Hard;
        }
    }
    ensure_equals(l1Hard, 0);
    ensure_equals(l2Hard, 0);
}

// P1: closed-form length must beat (or stay within 15% of) toLinear densify.
template<>
template<>
void object::test<9>()
{
    set_test_name("P1 length at most 1.15x toLinear");
    std::mt19937 rng(2817130497u ^ 81u);
    std::uniform_real_distribution<double> box(-100.0, 100.0);
    constexpr int nSample = 4000;
    std::vector<ExactCircularArc> sample;
    sample.reserve(static_cast<std::size_t>(nSample));
    while (static_cast<int>(sample.size()) < nSample) {
        ExactCircularArc a(CoordinateXY{box(rng), box(rng)},
                           CoordinateXY{box(rng), box(rng)},
                           CoordinateXY{box(rng), box(rng)});
        if (a.isArc()) {
            sample.push_back(a);
        }
    }
    for (int i = 0; i < 64; ++i) {
        (void)sample[static_cast<std::size_t>(i)].length();
        (void)sample[static_cast<std::size_t>(i)].toLinear(0.01);
    }
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    double sink = 0.0;
    for (const auto& a : sample) {
        sink += a.length();
    }
    const auto t1 = clock::now();
    for (const auto& a : sample) {
        sink += a.toLinear(0.01)->getLength();
    }
    const auto t2 = clock::now();
    const auto aNs = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    const auto dNs = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    const double ratio = dNs == 0 ? 0.0 : static_cast<double>(aNs) / static_cast<double>(dNs);
    std::cout << "P1 length/toLinear ns " << aNs << "/" << dNs
              << " = " << ratio << std::endl;
    ensure(sink != 0.0);
    ensure(ratio <= 1.15);
}

template<>
template<>
void object::test<10>()
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

} // namespace tut
