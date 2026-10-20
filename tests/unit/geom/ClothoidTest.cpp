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
#include <geos/io/ParseException.h>
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
using geos::geom::ClothoidReferenceLocation;
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

    static ClothoidReferenceLocation identityAt(double x, double y, double z = 0.0)
    {
        return ClothoidReferenceLocation{x, y, z, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
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
    set_test_name("WKT round trip of an arc-length pose uses the ISO body");

    auto curve = make(1.0, 2.0, 0.25, 0.0, 0.125, 8.0);
    geos::io::WKTWriter writer;
    writer.setRoundingPrecision(15);
    const std::string text = writer.write(curve.get());
    ensure(text.find("CLOTHOID ") == 0);
    ensure(text.find("REFERENCELOCATION") != std::string::npos);
    ensure(text.find("AFFINEPLACEMENT") != std::string::npos);
    ensure(text.find("SCALEFACTOR") != std::string::npos);
    ensure(text.find("STARTDISTANCE") != std::string::npos);
    ensure(text.find("ENDDISTANCE") != std::string::npos);
    ensure(text.find("CLOTHOID (1") == std::string::npos);

    auto back = reader_.read<Clothoid>(text);
    ensure(curve->equalsExact(back.get(), 1e-9));
    ensure_equals("value", back->getStartTangent(), curve->getStartTangent(), 1e-9);
    ensure_equals("value", back->getStartCurvature(), 0.0, 1e-9);
    ensure_equals("value", back->getEndCurvature(), 0.125, 1e-9);
    ensure_equals("value", back->getLength(), curve->getLength(), 1e-9);
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
    auto curve = reader_.read<Clothoid>(
        "CLOTHOID (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
        "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
        "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25)");
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
    ensure_equals("value", restored->getStartTangent(), curve->getStartTangent(), 1e-12);
    ensure_equals("value", restored->getScaleFactor(), curve->getScaleFactor(), 0.0);
    ensure_equals("value", restored->getStartDistance(), curve->getStartDistance(), 0.0);
    ensure_equals("value", restored->getEndDistance(), curve->getEndDistance(), 0.0);
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
    set_test_name("constant curvature is not an ST_Clothoid");

    // ISO/IEC 13249-3:2016 §4.2.11 places a Fresnel clothoid with a finite
    // scale factor. A straight chord (k0 = k1 = 0) has no such factor.
    // Clothoid DE-9IM is still not implemented.
    ensure_THROW(make(0.0, 0.0, 0.0, 0.0, 0.0, 2.0), geos::util::IllegalArgumentException);

    auto gentle = factory_->createClothoid(identityAt(0.0, 0.0), 1.0e6, 0.0, 2.0, false, false, 0.0, 0.0);
    ensure(gentle->hasCurvedComponents());
    ensure_equals("value", gentle->getLength(), 2.0, 1e-9);
    ensure_equals("value", gentle->getEndCoordinate().x, 2.0, 1e-6);
    ensure(std::abs(gentle->getEndCoordinate().y) < 1e-6);
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

template<>
template<>
void object::test<9>()
{
    set_test_name("ISO WKT and WKB round trip is bitwise for a short literal");

    const char* wkt =
        "CLOTHOID (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
        "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
        "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25)";
    auto curve = reader_.read<Clothoid>(wkt);
    ensure_equals("value", curve->getScaleFactor(), 50.0, 0.0);
    ensure_equals("value", curve->getStartDistance(), 0.0, 0.0);
    ensure_equals("value", curve->getEndDistance(), 25.0, 0.0);
    ensure_equals("value", curve->getReferenceLocation().x, 0.0, 0.0);
    ensure_equals("value", curve->getReferenceLocation().ux, 1.0, 0.0);
    ensure_equals("value", curve->getReferenceLocation().vy, 1.0, 0.0);

    geos::io::WKTWriter writer;
    const std::string text = writer.write(curve.get());
    auto back = reader_.read<Clothoid>(text);
    ensure(curve->equalsIdentical(back.get()));

    const char* hexLiteral =
        "0116000000"
        "0166000000"
        "0000000000000000"
        "0000000000000000"
        "02000000"
        "000000000000F03F"
        "0000000000000000"
        "0000000000000000"
        "000000000000F03F"
        "0000000000004940"
        "0000000000000000"
        "0000000000003940";
    geos::io::WKBWriter wkbWriter;
    std::stringstream hex;
    wkbWriter.writeHEX(*curve, hex);
    ensure_equals(hex.str(), std::string(hexLiteral));

    geos::io::WKBReader wkbReader(*factory_);
    std::stringstream hexIn(hexLiteral);
    auto fromHex = wkbReader.readHEX(hexIn);
    ensure(curve->equalsIdentical(fromHex.get()));
}

template<>
template<>
void object::test<10>()
{
    set_test_name("old arc-length CLOTHOID body and WKB type 18 are rejected");

    try {
        (void) reader_.read("CLOTHOID (0 0, 0, 0, 0.01, 25)");
        fail("expected ParseException");
    }
    catch (const geos::io::ParseException& ex) {
        const std::string msg(ex.what());
        ensure(msg.find("5.1.67") != std::string::npos);
        ensure(msg.find("not accepted") != std::string::npos);
    }

    geos::io::WKBReader wkbReader(*factory_);
    auto expectCircle = [&](const char* hexText) {
        try {
            std::stringstream hex(hexText);
            (void) wkbReader.readHEX(hex);
            fail(hexText);
        }
        catch (const geos::io::ParseException& ex) {
            const std::string msg(ex.what());
            ensure(msg.find("ST_Circle") != std::string::npos);
            ensure(msg.find("22") != std::string::npos);
        }
    };
    expectCircle("0112000000");
    expectCircle("01FA030000");

    try {
        std::stringstream hex(
            "0116000000"
            "0000000000000000"
            "0000000000000000"
            "0000000000000000"
            "0000000000000000"
            "0000000000000000"
            "0000000000000000");
        (void) wkbReader.readHEX(hex);
        fail("expected rejection of an arc-length body under type 22");
    }
    catch (const geos::io::ParseException& ex) {
        const std::string msg(ex.what());
        ensure(msg.find("5.1.68") != std::string::npos);
    }

    try {
        std::stringstream hex("0117000000");
        (void) wkbReader.readHEX(hex);
        fail("expected SPIRALCURVE refusal");
    }
    catch (const geos::io::ParseException& ex) {
        const std::string msg(ex.what());
        ensure(msg.find("not optional") != std::string::npos);
    }
}

template<>
template<>
void object::test<11>()
{
    set_test_name("Z, M, and ZM ISO forms round trip, and EMPTY is a header");

    const char* zWkt =
        "CLOTHOID Z (REFERENCELOCATION AFFINEPLACEMENT Z (LOCATION Z (1 2 9), "
        "REFERENCEDIRECTIONS (VECTOR Z (1 0 0), VECTOR Z (0 1 0))), "
        "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25)";
    auto zCurve = reader_.read<Clothoid>(zWkt);
    ensure(zCurve->hasZ());
    ensure(!zCurve->hasM());
    ensure_equals("value", zCurve->getCoordinatesRO()->getZ(0), 9.0, 0.0);
    ensure_equals("value", zCurve->getCoordinatesRO()->getZ(1), 9.0, 0.0);
    geos::io::WKTWriter writer;
    auto zBack = reader_.read<Clothoid>(writer.write(zCurve.get()));
    ensure(zCurve->equalsIdentical(zBack.get()));

    const char* mWkt =
        "CLOTHOID M (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
        "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
        "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25, STARTM 4, ENDM 9)";
    auto mCurve = reader_.read<Clothoid>(mWkt);
    ensure(mCurve->hasM());
    ensure(!mCurve->hasZ());
    ensure_equals("value", mCurve->getStartM(), 4.0, 0.0);
    ensure_equals("value", mCurve->getEndM(), 9.0, 0.0);
    ensure_equals("value", mCurve->getCoordinatesRO()->getM(0), 4.0, 0.0);
    ensure_equals("value", mCurve->getCoordinatesRO()->getM(1), 9.0, 0.0);
    auto mBack = reader_.read<Clothoid>(writer.write(mCurve.get()));
    ensure(mCurve->equalsIdentical(mBack.get()));

    auto line = mCurve->getLinearized(CurveToLineParams::maxDeviation(1.0));
    ensure_equals("value", line->getCoordinatesRO()->getM(0), 4.0, 1e-12);
    ensure_equals("value", line->getCoordinatesRO()->getM(line->getNumPoints() - 1), 9.0, 1e-12);

    const char* zmWkt =
        "CLOTHOID ZM (REFERENCELOCATION AFFINEPLACEMENT Z (LOCATION Z (1 2 3), "
        "REFERENCEDIRECTIONS (VECTOR Z (1 0 0), VECTOR Z (0 1 0))), "
        "SCALEFACTOR 8, STARTDISTANCE -2, ENDDISTANCE 2, STARTM 1, ENDM 5)";
    auto zmCurve = reader_.read<Clothoid>(zmWkt);
    ensure(zmCurve->hasZ());
    ensure(zmCurve->hasM());
    geos::io::WKBWriter wkbWriter;
    std::stringstream zmHex;
    wkbWriter.writeHEX(*zmCurve, zmHex);
    geos::io::WKBReader wkbReader(*factory_);
    auto zmBack = wkbReader.readHEX(zmHex);
    ensure(zmCurve->equalsIdentical(zmBack.get()));

    auto empty = factory_->createClothoid(false, false);
    std::stringstream emptyHex;
    wkbWriter.writeHEX(*empty, emptyHex);
    ensure_equals(emptyHex.str(), std::string("0116000000"));
    std::stringstream emptyZ;
    wkbWriter.writeHEX(*factory_->createClothoid(true, false), emptyZ);
    ensure_equals(emptyZ.str(), std::string("01FE030000"));
    std::stringstream emptyM;
    wkbWriter.writeHEX(*factory_->createClothoid(false, true), emptyM);
    ensure_equals(emptyM.str(), std::string("01E6070000"));
    std::stringstream emptyZM;
    wkbWriter.writeHEX(*factory_->createClothoid(true, true), emptyZM);
    ensure_equals(emptyZM.str(), std::string("01CE0B0000"));

    auto parsedMEmpty = reader_.read("CLOTHOID M EMPTY");
    ensure(parsedMEmpty->isEmpty());
    ensure(parsedMEmpty->hasM());
    ensure(!parsedMEmpty->hasZ());

    try {
        (void) reader_.read(
            "CLOTHOID (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
            "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
            "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25, STARTM 1, ENDM 2)");
        fail("expected measures to require the M flag");
    }
    catch (const geos::io::ParseException& ex) {
        const std::string msg(ex.what());
        ensure(msg.find("STARTM") != std::string::npos);
    }
    try {
        (void) reader_.read(
            "CLOTHOID M (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
            "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
            "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25)");
        fail("expected M to require STARTM and ENDM");
    }
    catch (const geos::io::ParseException& ex) {
        const std::string msg(ex.what());
        ensure(msg.find("STARTM") != std::string::npos);
    }
}

template<>
template<>
void object::test<12>()
{
    set_test_name("sheared placement and a span through the inflection");

    ClothoidReferenceLocation sheared{0, 0, 0, 2, 0, 0, 0, 1, 0};
    auto curve = factory_->createClothoid(sheared, 1.0, 0.0, 1.0, false, false, 0.0, 0.0);
    Parameters local;
    local.x0 = 0.0;
    local.y0 = 0.0;
    local.theta0 = 0.0;
    local.k0 = 0.0;
    local.k1 = 1.0;
    local.length = 1.0;
    double x = 0.0;
    double y = 0.0;
    position(local, 1.0, x, y);
    ensure_equals("value", curve->getEndCoordinate().x, 2.0 * x, 1e-9);
    ensure_equals("value", curve->getEndCoordinate().y, y, 1e-9);

    auto through = factory_->createClothoid(identityAt(0.0, 0.0), 10.0, -5.0, 5.0, false, false, 0.0, 0.0);
    ensure_equals("value", through->getLength(), 10.0, 1e-12);
    ensure_equals("value", through->getStartCoordinate().x, -through->getEndCoordinate().x, 1e-9);
    ensure_equals("value", through->getStartCoordinate().y, -through->getEndCoordinate().y, 1e-9);
    ensure(through->getEnvelopeInternal()->covers(0.0, 0.0));

    auto zero = factory_->createClothoid(identityAt(3.0, 4.0), 10.0, 1.0, 1.0, false, false, 0.0, 0.0);
    ensure(!zero->isEmpty());
    ensure_equals("value", zero->getLength(), 0.0, 0.0);
    ensure_equals("value", zero->getStartCoordinate().x, zero->getEndCoordinate().x, 0.0);
}

} // namespace tut
