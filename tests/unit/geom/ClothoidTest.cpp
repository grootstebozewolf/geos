#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>

#include <geos/algorithm/Clothoid.h>
#include <geos/algorithm/CurveToLineParams.h>
#include <geos/constants.h>
#include <geos/geom/Clothoid.h>
#include <geos/geom/CompoundCurve.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/io/WKBReader.h>
#include <geos/io/WKBWriter.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/util/IllegalArgumentException.h>
#include <geos/util/UnsupportedOperationException.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using geos::algorithm::CurveToLineParams;
using geos::algorithm::clothoid::Parameters;
using geos::algorithm::clothoid::position;
using geos::geom::Clothoid;
using geos::geom::CompoundCurve;
using geos::geom::CoordinateXY;
using geos::geom::CoordinateXYZM;
using geos::geom::SimpleCurve;

namespace tut {

struct test_clothoid_data {
    geos::geom::GeometryFactory::Ptr factory_ = geos::geom::GeometryFactory::create();
    geos::io::WKTReader reader_;

    std::unique_ptr<Clothoid> make(
        double x, double y, double theta0, double k0, double k1, double length) const
    {
        return factory_->createClothoid(
            CoordinateXYZM(x, y, geos::DoubleNotANumber, geos::DoubleNotANumber),
            false, false, theta0, k0, k1, length);
    }
};

typedef test_group<test_clothoid_data> group;
typedef group::object object;

group test_clothoid_group("geos::geom::Clothoid");

template<>
template<>
void object::test<1>()
{
    set_test_name("empty clothoid WKT and length");

    auto empty = factory_->createClothoid(false, false);
    ensure(empty->isEmpty());
    ensure_equals(empty->getLength(), 0.0);
    ensure_equals(empty->getNumPoints(), 0u);
    ensure(empty->hasCurvedTypes());
    ensure(!empty->hasCurvedComponents());
    ensure(empty->getEnvelopeInternal()->isNull());

    auto parsed = reader_.read("CLOTHOID EMPTY");
    ensure_equals(parsed->getGeometryType(), "Clothoid");
    ensure_equals(parsed->getGeometryTypeId(), geos::geom::GEOS_CLOTHOID);
    ensure(parsed->isEmpty());

    auto zempty = reader_.read("CLOTHOID Z EMPTY");
    ensure(zempty->isEmpty());
    ensure(zempty->hasZ());
    ensure(!zempty->hasM());
}

template<>
template<>
void object::test<2>()
{
    set_test_name("WKT round trip keeps the stored parameters");

    auto curve = make(1.0, 2.0, 0.25, 0.0, 0.125, 8.0);
    geos::io::WKTWriter writer;
    writer.setRoundingPrecision(15);
    const std::string text = writer.write(curve.get());
    ensure(text.find("CLOTHOID") == 0);

    auto back = reader_.read<Clothoid>(text);
    ensure(curve->equalsIdentical(back.get()));
    ensure_equals("value", back->getStartTangent(), 0.25, 0.0);
    ensure_equals("value", back->getStartCurvature(), 0.0, 0.0);
    ensure_equals("value", back->getEndCurvature(), 0.125, 0.0);
    ensure_equals("value", back->getLength(), 8.0, 0.0);
}

template<>
template<>
void object::test<3>()
{
    set_test_name("WKB round trip of a Z clothoid and of EMPTY");

    auto curve = factory_->createClothoid(
        CoordinateXYZM(1.5, -2.5, 9.0, geos::DoubleNotANumber),
        true, false, 0.3, -0.01, 0.02, 4.0);
    ensure(curve->hasZ());

    geos::io::WKBWriter writer;
    std::stringstream hex;
    writer.writeHEX(*curve, hex);
    geos::io::WKBReader wkbReader(*factory_);
    auto back = wkbReader.readHEX(hex);
    ensure_equals(back->getGeometryTypeId(), geos::geom::GEOS_CLOTHOID);
    ensure(curve->equalsIdentical(back.get()));
    ensure(back->hasZ());
    const auto* backClothoid = dynamic_cast<const Clothoid*>(back.get());
    ensure(backClothoid != nullptr);
    ensure_equals("value", backClothoid->getCoordinatesRO()->getZ(0), 9.0, 0.0);
    ensure_equals("value", backClothoid->getCoordinatesRO()->getZ(1), 9.0, 0.0);

    auto empty = factory_->createClothoid(false, false);
    std::stringstream emptyHex;
    writer.writeHEX(*empty, emptyHex);
    auto emptyBack = wkbReader.readHEX(emptyHex);
    ensure(emptyBack->isEmpty());
    ensure_equals(emptyBack->getGeometryTypeId(), geos::geom::GEOS_CLOTHOID);
}

template<>
template<>
void object::test<4>()
{
    set_test_name("IFC clothoid end point, length, and envelope");

    // tools/cad/fixtures/ifc_clothoid.json: k0 = 0, end radius 100, length 25.
    auto curve = reader_.read<Clothoid>("CLOTHOID (0 0, 0, 0, 0.01, 25)");
    ensure_equals("value", curve->getLength(), 25.0, 0.0);
    ensure_equals("value", curve->getEndCoordinate().x, 24.960965746830126, 1e-9);
    ensure_equals("value", curve->getEndCoordinate().y, 1.0405046700886817, 1e-9);
    ensure_equals("value", curve->getEndTangent(), 0.125, 1e-12);

    // Heading passes pi/2, so the box is wider than the endpoint chord.
    auto turning = make(0.0, 0.0, 0.0, 0.0, 2.0, geos::MATH_PI);
    const double sStar = geos::MATH_PI / std::sqrt(2.0);
    double xStar = 0.0;
    double yStar = 0.0;
    Parameters params;
    params.x0 = 0.0;
    params.y0 = 0.0;
    params.theta0 = 0.0;
    params.k0 = 0.0;
    params.k1 = 2.0;
    params.length = geos::MATH_PI;
    position(params, sStar, xStar, yStar);

    const geos::geom::Envelope* env = turning->getEnvelopeInternal();
    ensure(env->covers(xStar, yStar));
    ensure(env->getMaxX() > turning->getEndCoordinate().x + 0.05);
    ensure(env->covers(&turning->getStartCoordinate()));
    ensure(env->covers(&turning->getEndCoordinate()));
}

template<>
template<>
void object::test<5>()
{
    set_test_name("linearize and reverse agree on sampled positions");

    auto curve = make(0.0, 0.0, 0.2, 0.1, -0.05, 3.0);
    auto line = curve->getLinearized(CurveToLineParams::maxDeviation(1e-4));
    ensure_equals(line->getGeometryTypeId(), geos::geom::GEOS_LINESTRING);
    ensure(line->getNumPoints() >= 3u);
    ensure_equals("value", line->getCoordinatesRO()->front<CoordinateXY>().x, curve->getStartCoordinate().x, 1e-12);
    ensure_equals("value", line->getCoordinatesRO()->back<CoordinateXY>().x, curve->getEndCoordinate().x, 1e-9);
    ensure_equals("value", line->getCoordinatesRO()->back<CoordinateXY>().y, curve->getEndCoordinate().y, 1e-9);

    auto reversed = curve->reverse();
    Parameters forward;
    forward.x0 = curve->getStartCoordinate().x;
    forward.y0 = curve->getStartCoordinate().y;
    forward.theta0 = curve->getStartTangent();
    forward.k0 = curve->getStartCurvature();
    forward.k1 = curve->getEndCurvature();
    forward.length = curve->getLength();
    Parameters backward;
    backward.x0 = reversed->getStartCoordinate().x;
    backward.y0 = reversed->getStartCoordinate().y;
    backward.theta0 = reversed->getStartTangent();
    backward.k0 = reversed->getStartCurvature();
    backward.k1 = reversed->getEndCurvature();
    backward.length = reversed->getLength();

    const double s = 1.25;
    double x = 0.0;
    double y = 0.0;
    double xr = 0.0;
    double yr = 0.0;
    position(forward, s, x, y);
    position(backward, forward.length - s, xr, yr);
    ensure_equals("value", x, xr, 1e-8);
    ensure_equals("value", y, yr, 1e-8);

    auto restored = reversed->reverse();
    ensure_equals("value", restored->getStartTangent(), curve->getStartTangent() + 2.0 * geos::MATH_PI, 1e-9);
    double xe = 0.0;
    double ye = 0.0;
    Parameters again = forward;
    again.theta0 = restored->getStartTangent();
    again.x0 = restored->getStartCoordinate().x;
    again.y0 = restored->getStartCoordinate().y;
    position(again, s, xe, ye);
    ensure_equals("value", x, xe, 1e-8);
    ensure_equals("value", y, ye, 1e-8);
}

template<>
template<>
void object::test<6>()
{
    set_test_name("straight chord length and envelope, without a relate matrix");

    // oracle/de9im_clothoid_vectors.txt uses the chord (0, 0)-(2, 0).
    // This checks that degenerate curvature is that chord. The matrix in
    // the oracle is not evaluated; clothoid DE-9IM is not implemented.
    auto chord = make(0.0, 0.0, 0.0, 0.0, 0.0, 2.0);
    ensure(!chord->hasCurvedComponents());
    ensure(chord->hasCurvedTypes());
    ensure_equals("value", chord->getLength(), 2.0, 0.0);
    ensure_equals("value", chord->getEndCoordinate().x, 2.0, 1e-12);
    ensure_equals("value", chord->getEndCoordinate().y, 0.0, 1e-12);
    const geos::geom::Envelope* env = chord->getEnvelopeInternal();
    ensure_equals("value", env->getMinX(), 0.0, 0.0);
    ensure_equals("value", env->getMaxX(), 2.0, 1e-12);
    ensure_equals("value", env->getMinY(), 0.0, 1e-12);
    ensure_equals("value", env->getMaxY(), 0.0, 1e-12);
}

template<>
template<>
void object::test<7>()
{
    set_test_name("Halley chord placement and curved convexHull refusal");

    auto placed = factory_->createClothoidFromChord(
        CoordinateXY(166473.5960999988, 577235.5036999993),
        CoordinateXY(166558.2337500006, 577204.9115500003),
        -0.0006369426751592356, 0.0);
    ensure_equals("value", placed->getLength(), 90.0000162557141, 1e-6);
    ensure_equals("value", placed->getEndCoordinate().x, 166558.2337500006, 1e-4);
    ensure_equals("value", placed->getEndCoordinate().y, 577204.9115500003, 1e-4);

    auto curve = make(0.0, 0.0, 0.0, 0.0, 0.01, 25.0);
    try {
        (void) curve->convexHull();
        fail("expected UnsupportedOperationException");
    }
    catch (const geos::util::UnsupportedOperationException& ex) {
        const std::string msg(ex.what());
        ensure(msg.find("getLinearized") != std::string::npos);
    }

    ensure_THROW(
        factory_->createClothoid(
            CoordinateXYZM(0.0, 0.0, geos::DoubleNotANumber, geos::DoubleNotANumber),
            false, false, 0.0, 0.0, 0.0, -1.0),
        geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<8>()
{
    set_test_name("compound curve section linearizes");

    auto clothoid = make(0.0, 0.0, 0.0, 0.0, 0.01, 25.0);
    const CoordinateXY end = clothoid->getEndCoordinate();
    std::vector<std::unique_ptr<SimpleCurve>> sections;
    sections.push_back(std::move(clothoid));
    sections.push_back(factory_->createLineString({
        end,
        CoordinateXY(end.x, end.y + 1.0)
    }));
    std::unique_ptr<CompoundCurve> compound = factory_->createCompoundCurve(std::move(sections));
    auto line = compound->getLinearized(CurveToLineParams::maxDeviation(0.05));
    ensure_equals(line->getGeometryTypeId(), geos::geom::GEOS_LINESTRING);
    ensure(line->getNumPoints() > 3u);
    ensure_equals("value", line->getCoordinatesRO()->front<CoordinateXY>().x, 0.0, 1e-12);
    ensure_equals("value", line->getCoordinatesRO()->back<CoordinateXY>().y, end.y + 1.0, 1e-8);
}

} // namespace tut
