#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>

#include <geos/algorithm/CurveToLineParams.h>
#include <geos/algorithm/PointLocator.h>
#include <geos/constants.h>
#include <geos/geom/CircularString.h>
#include <geos/geom/Clothoid.h>
#include <geos/geom/CompoundCurve.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
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
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using geos::algorithm::CurveToLineParams;
using geos::algorithm::PointLocator;
using geos::geom::Clothoid;
using geos::geom::ClothoidReferenceLocation;
using geos::geom::CompoundCurve;
using geos::geom::Coordinate;
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
        return ClothoidReferenceLocation{
            Coordinate(x, y, z),
            Coordinate(1.0, 0.0, 0.0),
            Coordinate(0.0, 1.0, 0.0)};
    }

    std::unique_ptr<CompoundCurve> lineThenClothoid() const
    {
        std::vector<std::unique_ptr<SimpleCurve>> sections;
        sections.push_back(factory_->createLineString({
            CoordinateXY(0.0, 0.0), CoordinateXY(1.0, 0.0)}));
        sections.push_back(make(1.0, 0.0, 0.0, 0.0, 0.01, 5.0));
        return factory_->createCompoundCurve(std::move(sections));
    }
};

typedef test_group<test_clothoid_data> group;
typedef group::object object;

group test_clothoid_group("geos::geom::Clothoid");

template<>
template<>
void object::test<1>()
{
    set_test_name("type id");

    ensure_equals(static_cast<int>(geos::geom::GEOS_CLOTHOID), 13);
    auto empty = factory_->createClothoid(false, false);
    ensure_equals(empty->getGeometryTypeId(), geos::geom::GEOS_CLOTHOID);
    ensure_equals(empty->getGeometryType(), "Clothoid");
}

template<>
template<>
void object::test<2>()
{
    set_test_name("empty");

    auto empty = factory_->createClothoid(false, false);
    ensure(empty->isEmpty());
    ensure_equals(empty->getLength(), 0.0);
    ensure_equals(empty->getNumPoints(), 0u);
    ensure(empty->hasCurvedTypes());
    ensure(!empty->hasCurvedComponents());
    ensure(empty->getEnvelopeInternal()->isNull());

    auto parsed = reader_.read("CLOTHOID EMPTY");
    ensure(parsed->isEmpty());
    ensure_equals(parsed->getGeometryTypeId(), geos::geom::GEOS_CLOTHOID);

    auto zempty = reader_.read("CLOTHOID Z EMPTY");
    ensure(zempty->hasZ());
    ensure(!zempty->hasM());
    auto mempty = reader_.read("CLOTHOID M EMPTY");
    ensure(mempty->hasM());
    ensure(!mempty->hasZ());
}

template<>
template<>
void object::test<3>()
{
    set_test_name("WKT round trip");

    auto curve = make(1.0, 2.0, 0.25, 0.0, 0.125, 8.0);
    geos::io::WKTWriter writer;
    writer.setRoundingPrecision(15);
    const std::string text = writer.write(curve.get());
    ensure(text.find("REFERENCELOCATION") != std::string::npos);
    ensure(text.find("CLOTHOID (1") == std::string::npos);

    auto back = reader_.read<Clothoid>(text);
    ensure(curve->equalsExact(back.get(), 1e-9));
    ensure_equals("heading", back->getStartHeading(), curve->getStartHeading(), 1e-9);
    ensure_equals("k0", back->getStartCurvature(), 0.0, 1e-9);
    ensure_equals("k1", back->getEndCurvature(), 0.125, 1e-9);
    ensure_equals("length", back->getLength(), curve->getLength(), 1e-9);
}

template<>
template<>
void object::test<4>()
{
    set_test_name("WKB Z");

    auto curve = factory_->createClothoid(
                     CoordinateXYZM(1.5, -2.5, 9.0, geos::DoubleNotANumber),
                     true, false, 0.3, -0.01, 0.02, 4.0);
    geos::io::WKBWriter writer;
    std::stringstream hex;
    writer.writeHEX(*curve, hex);
    geos::io::WKBReader wkbReader(*factory_);
    auto back = wkbReader.readHEX(hex);
    ensure(curve->equalsIdentical(back.get()));
    const auto* backClothoid = dynamic_cast<const Clothoid*>(back.get());
    ensure(backClothoid != nullptr);
    ensure_equals("z0", backClothoid->getCoordinatesRO()->getZ(0), 9.0, 0.0);
    ensure_equals("z1", backClothoid->getCoordinatesRO()->getZ(1), 9.0, 0.0);
}

template<>
template<>
void object::test<5>()
{
    set_test_name("length");

    auto curve = reader_.read<Clothoid>(
                     "CLOTHOID (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
                     "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
                     "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25)");
    ensure_equals("length", curve->getLength(), 25.0, 0.0);
    ensure_equals("x", curve->getEndCoordinate().x, 24.960965746830126, 1e-9);
    ensure_equals("y", curve->getEndCoordinate().y, 1.0405046700886817, 1e-9);
    ensure_equals("heading", curve->getEndHeading(), 0.125, 1e-12);
}

template<>
template<>
void object::test<6>()
{
    set_test_name("envelope");

    auto turning = make(0.0, 0.0, 0.0, 0.0, 2.0, geos::MATH_PI);
    auto line = turning->getLinearized(CurveToLineParams::maxDeviation(1e-4));
    const geos::geom::Envelope* env = turning->getEnvelopeInternal();
    const auto* seq = line->getCoordinatesRO();
    for (std::size_t i = 0; i < seq->size(); ++i) {
        const CoordinateXY pt = seq->getAt<CoordinateXY>(i);
        ensure(env->covers(pt.x, pt.y));
    }
    ensure(env->getMaxX() > turning->getEndCoordinate().x + 0.05);
}

template<>
template<>
void object::test<7>()
{
    set_test_name("getLinearized()");

    auto curve = make(0.0, 0.0, 0.2, 0.1, -0.05, 3.0);
    auto line = curve->getLinearized(CurveToLineParams::maxDeviation(1e-4));
    ensure_equals(line->getGeometryTypeId(), geos::geom::GEOS_LINESTRING);
    ensure(line->getNumPoints() >= 3u);
    ensure_equals("x0", line->getCoordinatesRO()->front<CoordinateXY>().x, curve->getStartCoordinate().x, 1e-12);
    ensure_equals("x1", line->getCoordinatesRO()->back<CoordinateXY>().x, curve->getEndCoordinate().x, 1e-9);
    ensure_equals("y1", line->getCoordinatesRO()->back<CoordinateXY>().y, curve->getEndCoordinate().y, 1e-9);
}

template<>
template<>
void object::test<8>()
{
    set_test_name("getLinearized() Z/M");

    auto curve = reader_.read<Clothoid>(
                     "CLOTHOID ZM (REFERENCELOCATION AFFINEPLACEMENT Z (LOCATION Z (1 2 3), "
                     "REFERENCEDIRECTIONS (VECTOR Z (1 0 0), VECTOR Z (0 1 0))), "
                     "SCALEFACTOR 8, STARTDISTANCE -2, ENDDISTANCE 2, STARTM 1, ENDM 5)");
    auto line = curve->getLinearized(CurveToLineParams::maxDeviation(0.1));
    const auto* seq = line->getCoordinatesRO();
    ensure_equals("z0", seq->getZ(0), 3.0, 1e-12);
    ensure_equals("z1", seq->getZ(seq->size() - 1), 3.0, 1e-12);
    ensure_equals("m0", seq->getM(0), 1.0, 1e-12);
    ensure_equals("m1", seq->getM(seq->size() - 1), 5.0, 1e-12);
}

template<>
template<>
void object::test<9>()
{
    set_test_name("reverse");

    auto curve = make(0.0, 0.0, 0.2, 0.1, -0.05, 3.0);
    auto reversed = curve->reverse();
    ensure_equals("x", reversed->getStartCoordinate().x, curve->getEndCoordinate().x, 1e-9);
    ensure_equals("y", reversed->getStartCoordinate().y, curve->getEndCoordinate().y, 1e-9);
    ensure_equals("k0", reversed->getStartCurvature(), -curve->getEndCurvature(), 1e-12);
    ensure_equals("k1", reversed->getEndCurvature(), -curve->getStartCurvature(), 1e-12);
    ensure_equals("length", reversed->getLength(), curve->getLength(), 1e-9);

    auto restored = reversed->reverse();
    ensure_equals("heading", restored->getStartHeading(), curve->getStartHeading(), 1e-12);
    ensure_equals("scale", restored->getScaleFactor(), curve->getScaleFactor(), 0.0);
    ensure_equals("t0", restored->getStartDistance(), curve->getStartDistance(), 0.0);
    ensure_equals("t1", restored->getEndDistance(), curve->getEndDistance(), 0.0);
}

template<>
template<>
void object::test<10>()
{
    set_test_name("measures");

    auto curve = factory_->createClothoid(
                     CoordinateXYZM(0.0, 0.0, geos::DoubleNotANumber, 2.0),
                     false, true, 0.0, 0.0, 0.01, 10.0, 7.0);
    ensure(curve->hasM());
    ensure_equals("startM", curve->getStartM(), 2.0, 0.0);
    ensure_equals("endM", curve->getEndM(), 7.0, 0.0);
    ensure_equals("m0", curve->getCoordinatesRO()->getM(0), 2.0, 0.0);
    ensure_equals("m1", curve->getCoordinatesRO()->getM(1), 7.0, 0.0);
}

template<>
template<>
void object::test<11>()
{
    set_test_name("constant curvature");

    ensure_THROW(make(0.0, 0.0, 0.0, 0.0, 0.0, 2.0), geos::util::IllegalArgumentException);
    ensure_THROW(make(0.0, 0.0, 0.0, 0.01, 0.01, 2.0), geos::util::IllegalArgumentException);
    ensure_THROW(make(0.0, 0.0, 0.0, 0.0, 0.01, -1.0), geos::util::IllegalArgumentException);
    ensure_THROW(
        make(0.0, 0.0, 0.0, 0.0, 0.01, std::numeric_limits<double>::infinity()),
        geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<12>()
{
    set_test_name("gentle curve");

    auto gentle = factory_->createClothoid(identityAt(0.0, 0.0), 1.0e6, 0.0, 2.0, false, false, 0.0, 0.0);
    ensure(gentle->hasCurvedComponents());
    ensure_equals("length", gentle->getLength(), 2.0, 1e-9);
    ensure_equals("x", gentle->getEndCoordinate().x, 2.0, 1e-6);
    ensure(std::abs(gentle->getEndCoordinate().y) < 1e-6);
}

template<>
template<>
void object::test<13>()
{
    set_test_name("large turn");

    ensure_THROW(make(0.0, 0.0, 0.0, 0.0, 1e5, 20000.0), geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<14>()
{
    set_test_name("large curvature");

    ensure_THROW(make(0.0, 0.0, 0.0, 0.0, 1.0, 1e308), geos::util::IllegalArgumentException);
    ensure_THROW(make(0.0, 0.0, 0.0, 0.0, 1e308, 1.0), geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<15>()
{
    set_test_name("compound linearize");

    auto compound = lineThenClothoid();
    const CoordinateXY end = compound->getCurveN(1)->getEndCoordinate();
    auto line = compound->getLinearized(CurveToLineParams::maxDeviation(0.05));
    ensure(line->getNumPoints() > 3u);
    ensure_equals("x0", line->getCoordinatesRO()->front<CoordinateXY>().x, 0.0, 1e-12);
    ensure_equals("x1", line->getCoordinatesRO()->back<CoordinateXY>().x, end.x, 1e-8);
    ensure_equals("y1", line->getCoordinatesRO()->back<CoordinateXY>().y, end.y, 1e-8);
}

template<>
template<>
void object::test<16>()
{
    set_test_name("distance");

    auto compound = lineThenClothoid();
    auto point = factory_->createPoint(CoordinateXY(0.0, 10.0));
    ensure_THROW(point->distance(compound.get()), geos::util::UnsupportedOperationException);
}

template<>
template<>
void object::test<17>()
{
    set_test_name("intersection");

    auto compound = lineThenClothoid();
    auto line = factory_->createLineString({
        CoordinateXY(0.0, 0.0), CoordinateXY(1.0, 0.0)});
    ensure_THROW(compound->intersection(line.get()), geos::util::UnsupportedOperationException);
}

template<>
template<>
void object::test<18>()
{
    set_test_name("point location");

    auto compound = lineThenClothoid();
    const CoordinateXY end = compound->getCurveN(1)->getEndCoordinate();
    PointLocator locator;
    const geos::geom::Geometry* geom = compound.get();
    ensure_THROW(locator.locate(end, geom), geos::util::UnsupportedOperationException);
}

template<>
template<>
void object::test<19>()
{
    set_test_name("WKB hex");

    const char* wkt =
        "CLOTHOID (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
        "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
        "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25)";
    auto curve = reader_.read<Clothoid>(wkt);
    ensure_equals("scale", curve->getScaleFactor(), 50.0, 0.0);
    ensure_equals("ux", curve->getReferenceLocation().u.x, 1.0, 0.0);
    ensure_equals("vy", curve->getReferenceLocation().v.y, 1.0, 0.0);

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
void object::test<20>()
{
    set_test_name("WKB empty");

    geos::io::WKBWriter wkbWriter;
    std::stringstream emptyHex;
    wkbWriter.writeHEX(*factory_->createClothoid(false, false), emptyHex);
    ensure_equals(emptyHex.str(), std::string("0116000000"));
    std::stringstream emptyZ;
    wkbWriter.writeHEX(*factory_->createClothoid(true, false), emptyZ);
    ensure_equals(emptyZ.str(), std::string("0116000080"));
    std::stringstream emptyM;
    wkbWriter.writeHEX(*factory_->createClothoid(false, true), emptyM);
    ensure_equals(emptyM.str(), std::string("0116000040"));
    std::stringstream emptyZM;
    wkbWriter.writeHEX(*factory_->createClothoid(true, true), emptyZM);
    ensure_equals(emptyZM.str(), std::string("01160000C0"));

    geos::io::WKBWriter isoWriter(
        4, geos::io::ByteOrderValues::ENDIAN_LITTLE, false, geos::io::WKBConstants::wkbIso);
    std::stringstream isoZ;
    isoWriter.writeHEX(*factory_->createClothoid(true, false), isoZ);
    ensure_equals(isoZ.str(), std::string("01FE030000"));
    std::stringstream isoM;
    isoWriter.writeHEX(*factory_->createClothoid(false, true), isoM);
    ensure_equals(isoM.str(), std::string("01E6070000"));
    std::stringstream isoZM;
    isoWriter.writeHEX(*factory_->createClothoid(true, true), isoZM);
    ensure_equals(isoZM.str(), std::string("01CE0B0000"));
}

template<>
template<>
void object::test<21>()
{
    set_test_name("rejected WKT");

    ensure_THROW(reader_.read("CLOTHOID (0 0, 0, 0, 0.01, 25)"), geos::io::ParseException);
    ensure_THROW(reader_.read(
                     "CLOTHOID (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
                     "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
                     "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25, STARTM 1, ENDM 2)"),
                 geos::io::ParseException);
    ensure_THROW(reader_.read(
                     "CLOTHOID M (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
                     "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
                     "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25)"),
                 geos::io::ParseException);
}

template<>
template<>
void object::test<22>()
{
    set_test_name("WKB type 18");

    geos::io::WKBReader wkbReader(*factory_);
    std::stringstream circle("0112000000");
    ensure_THROW(wkbReader.readHEX(circle), geos::io::ParseException);
    std::stringstream circleZ("01FA030000");
    ensure_THROW(wkbReader.readHEX(circleZ), geos::io::ParseException);
}

template<>
template<>
void object::test<23>()
{
    set_test_name("WKB type 23");

    geos::io::WKBReader wkbReader(*factory_);
    std::stringstream spiral("0117000000");
    ensure_THROW(wkbReader.readHEX(spiral), geos::io::ParseException);
}

template<>
template<>
void object::test<24>()
{
    set_test_name("sheared placement");

    ClothoidReferenceLocation sheared{
        Coordinate(0.0, 0.0, 0.0),
        Coordinate(2.0, 0.0, 0.0),
        Coordinate(0.0, 1.0, 0.0)};
    auto curve = factory_->createClothoid(sheared, 1.0, 0.0, 1.0, false, false, 0.0, 0.0);
    ensure_equals("x", curve->getEndCoordinate().x, 1.9505753764006892, 1e-9);
    ensure_equals("y", curve->getEndCoordinate().y, 0.16371404737570058, 1e-9);
}

template<>
template<>
void object::test<25>()
{
    set_test_name("inflection");

    auto through = factory_->createClothoid(identityAt(0.0, 0.0), 10.0, -5.0, 5.0, false, false, 0.0, 0.0);
    ensure_equals("length", through->getLength(), 10.0, 1e-12);
    ensure_equals("x", through->getStartCoordinate().x, -through->getEndCoordinate().x, 1e-9);
    ensure_equals("y", through->getStartCoordinate().y, -through->getEndCoordinate().y, 1e-9);
    ensure(through->getEnvelopeInternal()->covers(0.0, 0.0));

    auto zero = factory_->createClothoid(identityAt(3.0, 4.0), 10.0, 1.0, 1.0, false, false, 0.0, 0.0);
    ensure_equals("length", zero->getLength(), 0.0, 0.0);
    ensure_equals("x", zero->getStartCoordinate().x, zero->getEndCoordinate().x, 0.0);
}

template<>
template<>
void object::test<26>()
{
    set_test_name("decreasing curvature");

    auto easing = make(0.0, 0.0, 0.0, 0.02, 0.01, 50.0);
    ensure_equals("k0", easing->getStartCurvature(), 0.02, 1e-12);
    ensure_equals("k1", easing->getEndCurvature(), 0.01, 1e-12);
    ensure_equals("x0", easing->getStartCoordinate().x, 0.0, 1e-12);
    ensure_equals("y0", easing->getStartCoordinate().y, 0.0, 1e-12);
    ensure_equals("heading", easing->getStartHeading(), 0.0, 1e-12);
    ensure(easing->getEndCoordinate().x > 40.0);
    ensure(easing->getEndCoordinate().y > 10.0);

    auto turned = make(3.0, 4.0, 1.0, -0.002, -0.03, 40.0);
    ensure_equals("x0", turned->getStartCoordinate().x, 3.0, 1e-12);
    ensure_equals("y0", turned->getStartCoordinate().y, 4.0, 1e-12);
    ensure_equals("k0", turned->getStartCurvature(), -0.002, 1e-12);
    ensure_equals("k1", turned->getEndCurvature(), -0.03, 1e-12);
}

template<>
template<>
void object::test<27>()
{
    set_test_name("Z WKT");

    auto curve = reader_.read<Clothoid>(
                     "CLOTHOID Z (REFERENCELOCATION AFFINEPLACEMENT Z (LOCATION Z (1 2 9), "
                     "REFERENCEDIRECTIONS (VECTOR Z (1 0 0), VECTOR Z (0 1 0))), "
                     "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25)");
    ensure(curve->hasZ());
    ensure(!curve->hasM());
    ensure_equals("z", curve->getCoordinatesRO()->getZ(0), 9.0, 0.0);
    geos::io::WKTWriter writer;
    auto back = reader_.read<Clothoid>(writer.write(curve.get()));
    ensure(curve->equalsIdentical(back.get()));
}

template<>
template<>
void object::test<28>()
{
    set_test_name("M WKT");

    auto curve = reader_.read<Clothoid>(
                     "CLOTHOID M (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
                     "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
                     "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25, STARTM 4, ENDM 9)");
    ensure_equals("startM", curve->getStartM(), 4.0, 0.0);
    ensure_equals("endM", curve->getEndM(), 9.0, 0.0);
    geos::io::WKTWriter writer;
    auto back = reader_.read<Clothoid>(writer.write(curve.get()));
    ensure(curve->equalsIdentical(back.get()));
}

template<>
template<>
void object::test<29>()
{
    set_test_name("bad placement");

    geos::io::WKBReader wkbReader(*factory_);
    std::stringstream hex(
        "0116000000"
        "0000000000000000"
        "0000000000000000"
        "0000000000000000"
        "0000000000000000"
        "0000000000000000"
        "0000000000000000");
    ensure_THROW(wkbReader.readHEX(hex), geos::io::ParseException);
}

template<>
template<>
void object::test<15>()
{
    set_test_name("compound curve short form");

    const char* wkt =
        "COMPOUNDCURVE ("
        "(116414.353 411964.758, 116410.740 411976.388), "
        "CLOTHOID (0, 0.005, 48), "
        "CIRCULARSTRING (116394.687 412021.591, 116284.527 412126.266, 116132.940 412123.450), "
        "CLOTHOID (0.005, 0, 42), "
        "(116095.653 412104.165, 115842.170 411961.603))";

    auto curve = reader_.read<CompoundCurve>(wkt);
    ensure_equals(curve->getNumCurves(), 5u);

    const auto* entry = dynamic_cast<const Clothoid*>(curve->getCurveN(1));
    const auto* exit = dynamic_cast<const Clothoid*>(curve->getCurveN(3));
    ensure(entry != nullptr);
    ensure(exit != nullptr);
    ensure_equals("value", entry->getStartCurvature(), 0.0, 1e-12);
    ensure_equals("value", entry->getEndCurvature(), 0.005, 1e-12);
    ensure_equals("value", entry->getLength(), 48.0, 1e-9);
    ensure_equals("value", entry->getStartCoordinate().x, 116410.740, 1e-9);
    ensure_equals("value", entry->getStartCoordinate().y, 411976.388, 1e-9);
    ensure_equals("value", exit->getStartCurvature(), 0.005, 1e-12);
    ensure_equals("value", exit->getEndCurvature(), 0.0, 1e-12);
    ensure_equals("value", exit->getLength(), 42.0, 1e-9);

    const CoordinateXY& arcStart = curve->getCurveN(2)->getStartCoordinate();
    const CoordinateXY& lineStart = curve->getCurveN(4)->getStartCoordinate();
    ensure_equals("value", arcStart.x, 116394.687, 0.0);
    ensure_equals("value", arcStart.y, 412021.591, 0.0);
    ensure_equals("value", lineStart.x, 116095.653, 0.0);
    ensure_equals("value", lineStart.y, 412104.165, 0.0);
    ensure_equals("entry drift", entry->getEndCoordinate().distance(arcStart), 0.0, 2e-3);
    ensure_equals("exit drift", exit->getEndCoordinate().distance(lineStart), 0.0, 2e-3);

    auto linearized = curve->getLinearized(CurveToLineParams());
    auto contains = [&linearized](const CoordinateXY& point) {
        for (std::size_t i = 0; i < linearized->getNumPoints(); ++i) {
            if (linearized->getCoordinateN(i).equals2D(point)) {
                return true;
            }
        }
        return false;
    };
    ensure(contains(arcStart));
    ensure(contains(lineStart));

    geos::io::WKTWriter writer;
    const std::string text = writer.write(curve.get());
    ensure(text.find("CLOTHOID (0, 0.005, 48)") != std::string::npos);
    ensure(text.find("CLOTHOID (0.005, 0, 42)") != std::string::npos);
    ensure(text.find("REFERENCELOCATION") == std::string::npos);

    auto roundTrip = reader_.read<CompoundCurve>(text);
    ensure_equals("value", roundTrip->getCurveN(1)->getLength(), 48.0, 1e-6);
    ensure_equals("value", roundTrip->getCurveN(3)->getLength(), 42.0, 1e-6);
    ensure_equals("value", roundTrip->getCurveN(2)->getStartCoordinate().x, arcStart.x, 0.0);
    ensure_equals("value", roundTrip->getCurveN(2)->getStartCoordinate().y, arcStart.y, 0.0);

    geos::io::WKBWriter wkbWriter;
    std::stringstream hex;
    wkbWriter.writeHEX(*curve, hex);
    geos::io::WKBReader wkbReader(*factory_);
    std::stringstream hexIn(hex.str());
    auto fromWkb = wkbReader.readHEX(hexIn);
    auto* wkbCurve = dynamic_cast<CompoundCurve*>(fromWkb.get());
    ensure(wkbCurve != nullptr);
    ensure_equals(wkbCurve->getCurveN(1)->getGeometryTypeId(), geos::geom::GEOS_CLOTHOID);
    ensure_equals("value", wkbCurve->getCurveN(2)->getStartCoordinate().x, arcStart.x, 0.0);
    ensure_equals("entry drift",
        wkbCurve->getCurveN(1)->getEndCoordinate().distance(arcStart),
        entry->getEndCoordinate().distance(arcStart), 0.0);
    const std::string wkbText = writer.write(wkbCurve);
    ensure(wkbText.find("CLOTHOID (0, 0.005, 48)") != std::string::npos);
}

template<>
template<>
void object::test<16>()
{
    set_test_name("short form reverse");

    auto parsed = reader_.read<CompoundCurve>(
        "COMPOUNDCURVE ("
        "(116414.353 411964.758, 116410.740 411976.388), "
        "CLOTHOID (0, 0.005, 48), "
        "CIRCULARSTRING (116394.687 412021.591, 116284.527 412126.266, 116132.940 412123.450), "
        "CLOTHOID (0.005, 0, 42), "
        "(116095.653 412104.165, 115842.170 411961.603))");
    auto reversed = parsed->reverse();
    const auto* exit = dynamic_cast<const Clothoid*>(reversed->getCurveN(1));
    const auto* entry = dynamic_cast<const Clothoid*>(reversed->getCurveN(3));
    ensure(exit != nullptr);
    ensure(entry != nullptr);
    ensure_equals("value", exit->getStartCurvature(), 0.0, 1e-12);
    ensure_equals("value", exit->getEndCurvature(), -0.005, 1e-12);
    ensure_equals("value", exit->getLength(), 42.0, 1e-9);
    ensure_equals("value", entry->getStartCurvature(), -0.005, 1e-12);
    ensure_equals("value", entry->getEndCurvature(), 0.0, 1e-12);
    ensure_equals("value", entry->getLength(), 48.0, 1e-9);

    auto clothoid = make(0.0, 0.0, 0.0, 0.0, 0.01, 20.0);
    const double theta = clothoid->getEndTangent();
    const CoordinateXY end = clothoid->getEndCoordinate();
    const double radius = 100.0;
    const double cx = end.x - std::sin(theta) * radius;
    const double cy = end.y + std::cos(theta) * radius;
    const double ang0 = std::atan2(end.y - cy, end.x - cx);

    auto lineCoords = std::make_unique<geos::geom::CoordinateSequence>(static_cast<std::size_t>(2), false, false);
    lineCoords->setAt(CoordinateXY(-5.0, 0.0), 0);
    lineCoords->setAt(CoordinateXY(0.0, 0.0), 1);
    auto arcCoords = std::make_unique<geos::geom::CoordinateSequence>(static_cast<std::size_t>(3), false, false);
    arcCoords->setAt(end, 0);
    arcCoords->setAt(CoordinateXY(
        cx + radius * std::cos(ang0 + 15.0 * geos::MATH_PI / 180.0),
        cy + radius * std::sin(ang0 + 15.0 * geos::MATH_PI / 180.0)), 1);
    arcCoords->setAt(CoordinateXY(
        cx + radius * std::cos(ang0 + 30.0 * geos::MATH_PI / 180.0),
        cy + radius * std::sin(ang0 + 30.0 * geos::MATH_PI / 180.0)), 2);

    std::vector<std::unique_ptr<SimpleCurve>> parts;
    std::unique_ptr<SimpleCurve> line = factory_->createLineString(std::move(lineCoords));
    std::unique_ptr<SimpleCurve> spiral = std::move(clothoid);
    std::unique_ptr<SimpleCurve> arc = factory_->createCircularString(std::move(arcCoords));
    parts.push_back(std::move(line));
    parts.push_back(std::move(spiral));
    parts.push_back(std::move(arc));
    auto exact = factory_->createCompoundCurve(std::move(parts));

    geos::io::WKTWriter writer;
    const std::string text = writer.write(exact.get());
    ensure(text.find("CLOTHOID (0, 0.01, 20)") != std::string::npos);

    auto exactReversed = exact->reverse();
    const std::string reversedText = writer.write(exactReversed.get());
    ensure(reversedText.find("CLOTHOID (-0.01, 0, 20)") != std::string::npos);
    auto back = reader_.read<CompoundCurve>(reversedText);
    const auto* member = dynamic_cast<const Clothoid*>(back->getCurveN(1));
    ensure(member != nullptr);
    ensure_equals("value", member->getStartCurvature(), -0.01, 1e-12);
    ensure_equals("value", member->getEndCurvature(), 0.0, 1e-12);
    ensure_equals("value", member->getLength(), 20.0, 1e-9);
}

template<>
template<>
void object::test<17>()
{
    set_test_name("short form rejected outside a continuation");

    ensure_THROW(reader_.read("CLOTHOID (0, 0.005, 48)"), geos::io::ParseException);
    ensure_THROW(reader_.read("COMPOUNDCURVE (CLOTHOID (0, 0.005, 48), (0 0, 1 0))"), geos::io::ParseException);
    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0.01, 0.01, 5))"), geos::util::IllegalArgumentException);
    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0, 0.01, 0))"), geos::util::IllegalArgumentException);
    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0, 0.01, -3))"), geos::util::IllegalArgumentException);
    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (NaN, 0.01, 5))"), geos::util::IllegalArgumentException);
    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0, Inf, 5))"), geos::util::IllegalArgumentException);
    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), (1.1 0, 2 0))"), geos::util::IllegalArgumentException);

    auto mismatch = reader_.read<CompoundCurve>("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0.02, 0, 10))");
    ensure_equals(mismatch->getNumCurves(), 2u);
    geos::io::WKTWriter writer;
    const std::string text = writer.write(mismatch.get());
    ensure(text.find("REFERENCELOCATION") != std::string::npos);
    ensure(text.find("CLOTHOID (0.02") == std::string::npos);

    auto lifted = reader_.read<CompoundCurve>("COMPOUNDCURVE Z ((0 0 5, 3 0 5), CLOTHOID (0, 0.01, 6))");
    const auto* clothoid = dynamic_cast<const Clothoid*>(lifted->getCurveN(1));
    ensure(clothoid != nullptr);
    ensure(clothoid->hasZ());
    ensure_equals("value", clothoid->getCoordinatesRO()->getZ(0), 5.0, 0.0);
    ensure_equals("value", clothoid->getCoordinatesRO()->getZ(1), 5.0, 0.0);
    const std::string liftedText = writer.write(lifted.get());
    ensure(liftedText.find("CLOTHOID Z (0, 0.01, 6)") != std::string::npos);
}

} // namespace tut
