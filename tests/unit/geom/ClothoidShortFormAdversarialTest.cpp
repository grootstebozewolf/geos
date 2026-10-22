#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>

#include <geos/algorithm/CurveToLineParams.h>
#include <geos/geom/Clothoid.h>
#include <geos/geom/CompoundCurve.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/io/ParseException.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/util/IllegalArgumentException.h>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <memory>
#include <string>
#include <vector>

using geos::algorithm::CurveToLineParams;
using geos::geom::Clothoid;
using geos::geom::CompoundCurve;
using geos::geom::CoordinateXY;
using geos::geom::CoordinateXYZM;
using geos::geom::SimpleCurve;

namespace tut {

struct test_clothoid_short_form_data {
    geos::geom::GeometryFactory::Ptr factory_ = geos::geom::GeometryFactory::create();
    geos::io::WKTReader reader_;

    static bool hasPoint(const geos::geom::LineString& line, const CoordinateXY& point)
    {
        for (std::size_t i = 0; i < line.getNumPoints(); ++i) {
            if (line.getCoordinateN(i).equals2D(point)) {
                return true;
            }
        }
        return false;
    }
};

typedef test_group<test_clothoid_short_form_data> group;
typedef group::object object;

group test_clothoid_short_form_group("geos::geom::ClothoidShortFormAdversarial");

template<>
template<>
void object::test<1>()
{
    set_test_name("ISO clothoid misses the previous end");

    ensure_THROW(reader_.read(
        "COMPOUNDCURVE ((0 0, 1 0), "
        "CLOTHOID (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (10 10), "
        "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
        "SCALEFACTOR 1, STARTDISTANCE 0, ENDDISTANCE 1))"),
        geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<2>()
{
    set_test_name("ISO clothoid misses the next start");

    ensure_THROW(reader_.read(
        "COMPOUNDCURVE ("
        "CLOTHOID (REFERENCELOCATION AFFINEPLACEMENT (LOCATION (0 0), "
        "REFERENCEDIRECTIONS (VECTOR (1 0), VECTOR (0 1))), "
        "SCALEFACTOR 50, STARTDISTANCE 0, ENDDISTANCE 25), "
        "(0 0, 1 0))"),
        geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<3>()
{
    set_test_name("start curvature after a line");

    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0.02, 0, 10))"),
        geos::io::ParseException);
    auto curve = reader_.read<CompoundCurve>("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0, 0.01, 10))");
    const auto* clothoid = dynamic_cast<const Clothoid*>(curve->getCurveN(1));
    ensure(clothoid != nullptr);
    ensure(clothoid->getStartCoordinate().x == 1.0);
    ensure(clothoid->getStartCoordinate().y == 0.0);
}

template<>
template<>
void object::test<4>()
{
    set_test_name("start curvature after a CCW arc");

    ensure_THROW(reader_.read(
        "COMPOUNDCURVE (CIRCULARSTRING (1 0, 0 1, -1 0), CLOTHOID (0, 0.01, 5))"),
        geos::io::ParseException);
    auto curve = reader_.read<CompoundCurve>(
        "COMPOUNDCURVE (CIRCULARSTRING (1 0, 0 1, -1 0), CLOTHOID (1, 0, 5))");
    const auto* clothoid = dynamic_cast<const Clothoid*>(curve->getCurveN(1));
    ensure(clothoid != nullptr);
    ensure_equals("start curvature", clothoid->getStartCurvature(), 1.0, 1e-12);
    ensure(clothoid->getStartCoordinate().x == -1.0);
    ensure(clothoid->getStartCoordinate().y == 0.0);
}

template<>
template<>
void object::test<5>()
{
    set_test_name("start curvature after a CW arc");

    ensure_THROW(reader_.read(
        "COMPOUNDCURVE (CIRCULARSTRING (1 0, 0 -1, -1 0), CLOTHOID (0, -0.01, 5))"),
        geos::io::ParseException);
    auto curve = reader_.read<CompoundCurve>(
        "COMPOUNDCURVE (CIRCULARSTRING (1 0, 0 -1, -1 0), CLOTHOID (-1, 0, 5))");
    const auto* clothoid = dynamic_cast<const Clothoid*>(curve->getCurveN(1));
    ensure(clothoid != nullptr);
    ensure_equals("start curvature", clothoid->getStartCurvature(), -1.0, 1e-12);
}

template<>
template<>
void object::test<6>()
{
    set_test_name("start curvature after a clothoid");

    ensure_THROW(reader_.read(
        "COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0, 0.01, 10), CLOTHOID (0, 0.02, 8))"),
        geos::io::ParseException);
    auto curve = reader_.read<CompoundCurve>(
        "COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0, 0.01, 10), CLOTHOID (0.01, 0, 8))");
    const auto* second = dynamic_cast<const Clothoid*>(curve->getCurveN(2));
    ensure(second != nullptr);
    ensure_equals("start curvature", second->getStartCurvature(), 0.01, 1e-12);
    ensure(second->getStartCoordinate() == curve->getCurveN(1)->getEndCoordinate());
}

template<>
template<>
void object::test<7>()
{
    set_test_name("ZM round trip");

    auto curve = reader_.read<CompoundCurve>(
        "COMPOUNDCURVE ZM ((0 0 4 7, 1 0 4 7), CLOTHOID (0, 0.01, 5))");
    const auto* clothoid = dynamic_cast<const Clothoid*>(curve->getCurveN(1));
    ensure(clothoid != nullptr);
    ensure_equals("z", clothoid->getCoordinatesRO()->getZ(0), 4.0, 0.0);
    ensure_equals("m", clothoid->getCoordinatesRO()->getM(0), 7.0, 0.0);
    geos::io::WKTWriter writer;
    const std::string text = writer.write(curve.get());
    ensure(text.find("CLOTHOID ZM (0, 0.01, 5)") != std::string::npos);
    auto roundTrip = reader_.read<CompoundCurve>(text);
    const auto* again = dynamic_cast<const Clothoid*>(roundTrip->getCurveN(1));
    ensure(again != nullptr);
    ensure_equals("z", again->getCoordinatesRO()->getZ(1), 4.0, 0.0);
    ensure_equals("m", again->getCoordinatesRO()->getM(1), 7.0, 0.0);
}

template<>
template<>
void object::test<8>()
{
    set_test_name("different ZM writes ISO");

    auto lineCoords = std::make_unique<geos::geom::CoordinateSequence>(static_cast<std::size_t>(2), true, true);
    lineCoords->setAt(CoordinateXYZM(0.0, 0.0, 1.0, 2.0), 0);
    lineCoords->setAt(CoordinateXYZM(1.0, 0.0, 1.0, 2.0), 1);
    auto spiral = factory_->createClothoid(
        CoordinateXYZM(1.0, 0.0, 5.0, 9.0), true, true, 0.0, 0.0, 0.01, 5.0, 9.0);
    std::vector<std::unique_ptr<SimpleCurve>> parts;
    parts.push_back(factory_->createLineString(std::move(lineCoords)));
    parts.push_back(std::move(spiral));
    auto curve = factory_->createCompoundCurve(std::move(parts));
    geos::io::WKTWriter writer;
    const std::string text = writer.write(curve.get());
    ensure(text.find("REFERENCELOCATION") != std::string::npos);
    ensure(text.find("CLOTHOID ZM (0") == std::string::npos);
}

template<>
template<>
void object::test<9>()
{
    set_test_name("huge turn");

    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0, 1e5, 20000))"),
        geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<10>()
{
    set_test_name("one metre gap");

    auto spiral = factory_->createClothoid(
        CoordinateXYZM(1.0, 0.0, geos::DoubleNotANumber, geos::DoubleNotANumber),
        false, false, 0.0, 0.0, 0.01, 20.0);
    const CoordinateXY end = spiral->getEndCoordinate();
    std::ostringstream wkt;
    wkt << std::setprecision(17);
    wkt << "COMPOUNDCURVE ((0 0, 1 0), CLOTHOID (0, 0.01, 20), ("
        << (end.x + 1.0) << " " << end.y << ", "
        << (end.x + 2.0) << " " << end.y << "))";
    auto curve = reader_.read<CompoundCurve>(wkt.str());
    const auto* clothoid = dynamic_cast<const Clothoid*>(curve->getCurveN(1));
    ensure(clothoid != nullptr);
    const CoordinateXY typed = curve->getCurveN(2)->getStartCoordinate();
    ensure_equals("drift", clothoid->getEndCoordinate().distance(typed), 1.0, 1e-9);
    ensure(typed.x != clothoid->getEndCoordinate().x || typed.y != clothoid->getEndCoordinate().y);
    auto line = curve->getLinearized(CurveToLineParams());
    ensure(hasPoint(*line, typed));
    ensure(hasPoint(*line, clothoid->getEndCoordinate()));
}

template<>
template<>
void object::test<11>()
{
    set_test_name("CLOTHOID Z on an XY compound");

    ensure_THROW(reader_.read("COMPOUNDCURVE ((0 0, 1 0), CLOTHOID Z (0, 0.01, 5))"),
        geos::io::ParseException);
}

template<>
template<>
void object::test<12>()
{
    set_test_name("start curvature off the previous member by 5e-7");

    geos::io::WKTWriter writer;
    for (double k0 : {5e-7, 5e-8}) {
        auto lineCoords = std::make_unique<geos::geom::CoordinateSequence>(static_cast<std::size_t>(2), false, false);
        lineCoords->setAt(CoordinateXY(-5.0, 0.0), 0);
        lineCoords->setAt(CoordinateXY(0.0, 0.0), 1);
        std::vector<std::unique_ptr<SimpleCurve>> parts;
        parts.push_back(factory_->createLineString(std::move(lineCoords)));
        parts.push_back(factory_->createClothoid(
            CoordinateXYZM(0.0, 0.0, geos::DoubleNotANumber, geos::DoubleNotANumber),
            false, false, 0.0, k0, 0.01, 6.0));
        auto curve = factory_->createCompoundCurve(std::move(parts));

        const std::string text = writer.write(curve.get());
        const bool shortForm = text.find("REFERENCELOCATION") == std::string::npos;
        ensure_equals(text, shortForm, k0 < 1e-7);
        if (shortForm) {
            auto back = reader_.read<CompoundCurve>(text);
            const auto* clothoid = dynamic_cast<const Clothoid*>(back->getCurveN(1));
            ensure(clothoid != nullptr);
            ensure_equals("start curvature", clothoid->getStartCurvature(), k0, 1e-12);
        }
    }
}

template<>
template<>
void object::test<13>()
{
    set_test_name("reversed ProRail 823_12V_4.3 is written short");

    auto curve = reader_.read<CompoundCurve>(
        "COMPOUNDCURVE ("
        "(116414.353 411964.758, 116410.740 411976.388), "
        "CLOTHOID (0, 0.005, 48), "
        "CIRCULARSTRING (116394.687 412021.591, 116284.527 412126.266, 116132.940 412123.450), "
        "CLOTHOID (0.005, 0, 42), "
        "(116095.653 412104.165, 115842.170 411961.603))");
    auto reversed = curve->reverse();

    geos::io::WKTWriter writer;
    const std::string text = writer.write(reversed.get());
    ensure(text, text.find("CLOTHOID (0, -0.005, 42)") != std::string::npos);
    ensure(text, text.find("CLOTHOID (-0.005, 0, 48)") != std::string::npos);

    auto back = reader_.read<CompoundCurve>(text);
    ensure_equals("members", back->getNumCurves(), 5u);
    ensure_equals("start", back->getCurveN(1)->getStartCoordinate().x,
        reversed->getCurveN(0)->getEndCoordinate().x, 0.0);
    ensure_equals("drift", back->getCurveN(1)->getEndCoordinate().distance(
        reversed->getCurveN(1)->getEndCoordinate()), 0.0, 1e-3);
}

} // namespace tut
