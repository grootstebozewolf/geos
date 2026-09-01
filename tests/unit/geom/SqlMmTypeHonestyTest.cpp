// SQL/MM ISO/IEC 13249-3 type honesty (Proofs #660 leftovers).
// Named getLinearized is the only allowed linear fallback.
// Silent LINESTRING / POLYGON / MULTILINESTRING / MULTIPOLYGON emit is a refuse.

#include <tut/tut.hpp>
#include <utility.h>

#include <geos/algorithm/CurveToLineParams.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Surface.h>
#include <geos/io/WKTReader.h>
#include <geos/operation/overlayng/OverlayNG.h>
#include <geos/operation/union/CascadedPolygonUnion.h>
#include <geos/operation/union/UnionStrategy.h>
#include <geos/util/UnsupportedOperationException.h>

#include <memory>
#include <string>
#include <vector>

using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::geom::GeometryTypeId;
using geos::geom::Surface;
using geos::io::WKTReader;
using geos::operation::geounion::CascadedPolygonUnion;
using geos::operation::geounion::UnionStrategy;
using geos::operation::overlayng::OverlayNG;

namespace tut {

struct test_sqlmmtypehonesty_data {
    WKTReader r;

    std::unique_ptr<Geometry> read(const std::string& wkt)
    {
        return r.read(wkt);
    }
};

typedef test_group<test_sqlmmtypehonesty_data> group;
typedef group::object object;

group test_sqlmmtypehonesty_group("geos::geom::SqlMmTypeHonesty");

namespace {

class MixedDimLinearCurvePolygonUnion : public UnionStrategy {
public:
    std::unique_ptr<Geometry> Union(const Geometry*, const Geometry*) override
    {
        WKTReader reader;
        // All-linear CurvePolygon (no arcs) plus a leftover line.
        // restrictToSurfaces must not pack the CurvePolygon into MultiPolygon.
        return reader.read(
            "GEOMETRYCOLLECTION ("
            "CURVEPOLYGON (LINESTRING (0 0, 2 0, 2 2, 0 2, 0 0)), "
            "LINESTRING (10 10, 11 11))");
    }

    bool isFloatingPrecision() const override
    {
        return true;
    }
};

bool isLinearCollectionType(GeometryTypeId id)
{
    return id == GeometryTypeId::GEOS_LINESTRING
        || id == GeometryTypeId::GEOS_LINEARRING
        || id == GeometryTypeId::GEOS_POLYGON
        || id == GeometryTypeId::GEOS_MULTILINESTRING
        || id == GeometryTypeId::GEOS_MULTIPOLYGON;
}

} // namespace

// Ticket 36: all-linear CurvePolygon is still a curved type.
template<>
template<>
void object::test<1>()
{
    set_test_name("restrictToSurfaces keeps all-linear CurvePolygon as MultiSurface");

    MixedDimLinearCurvePolygonUnion strategy;
    auto dummy0 = read("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))");
    auto dummy1 = read("POLYGON ((2 2, 3 2, 3 3, 2 3, 2 2))");
    std::vector<const Surface*> polys;
    polys.push_back(static_cast<const Surface*>(dummy0.get()));
    polys.push_back(static_cast<const Surface*>(dummy1.get()));

    auto result = CascadedPolygonUnion::Union(polys, &strategy, nullptr);
    ensure(result != nullptr);
    ensure_equals(result->getGeometryTypeId(), GeometryTypeId::GEOS_MULTISURFACE);
    ensure_equals(result->getNumGeometries(), 1u);
    ensure_equals(result->getGeometryN(0)->getGeometryTypeId(), GeometryTypeId::GEOS_CURVEPOLYGON);
    ensure(result->getGeometryN(0)->hasCurvedTypes());
    ensure(!result->getGeometryN(0)->hasCurvedComponents());
}

// Ticket 30: hull refuses CompoundCurve and names getLinearized.
template<>
template<>
void object::test<2>()
{
    set_test_name("CompoundCurve convexHull refuses and names getLinearized");

    auto cc = read("COMPOUNDCURVE ((0 0, 1 0), CIRCULARSTRING (1 0, 1 1, 0 1))");
    try {
        (void) cc->convexHull();
        fail("expected UnsupportedOperationException");
    }
    catch (const geos::util::UnsupportedOperationException& ex) {
        const std::string msg(ex.what());
        ensure(msg.find("getLinearized") != std::string::npos);
    }
}

// Ticket 30: named getLinearized is the stamp (COMPOUNDCURVE → LINESTRING).
template<>
template<>
void object::test<3>()
{
    set_test_name("getLinearized is the named CompoundCurve stamp");

    auto cc = read("COMPOUNDCURVE ((0 0, 1 0), CIRCULARSTRING (1 0, 1 1, 0 1))");
    auto lin = cc->getLinearized(geos::algorithm::CurveToLineParams::stepSizeDegrees(45));
    ensure_equals(lin->getGeometryTypeId(), GeometryTypeId::GEOS_LINESTRING);
}

// Ticket 32: all-linear CurvePolygon overlay stays CurvePolygon, not Polygon.
template<>
template<>
void object::test<4>()
{
    set_test_name("all-linear CurvePolygon overlay does not emit Polygon");

    auto a = read("CURVEPOLYGON (LINESTRING (0 0, 2 0, 2 2, 0 2, 0 0))");
    auto b = read("CURVEPOLYGON (LINESTRING (1 1, 3 1, 3 3, 1 3, 1 1))");
    auto result = OverlayNG::overlay(a.get(), b.get(), OverlayNG::UNION);
    ensure(result != nullptr);
    ensure(result->hasCurvedTypes());
    ensure(!isLinearCollectionType(result->getGeometryTypeId()));
}

// Ticket 34: MultiCurve overlay stays a curve collection, not MultiLineString.
template<>
template<>
void object::test<5>()
{
    set_test_name("MultiCurve overlay does not emit MultiLineString");

    auto a = read("MULTICURVE (CIRCULARSTRING (0 0, 1 1, 2 0), (3 0, 4 0))");
    auto b = read("POINT (10 10)");
    auto result = OverlayNG::overlay(a.get(), b.get(), OverlayNG::UNION);
    ensure(result != nullptr);
    ensure(!isLinearCollectionType(result->getGeometryTypeId()));
    ensure(result->hasCurvedTypes());
}

// Ticket 36: MultiSurface overlay stays MultiSurface / CurvePolygon, not MultiPolygon.
template<>
template<>
void object::test<6>()
{
    set_test_name("MultiSurface overlay does not emit MultiPolygon");

    auto a = read("MULTISURFACE (CURVEPOLYGON (LINESTRING (0 0, 1 0, 1 1, 0 1, 0 0)))");
    auto b = read("MULTISURFACE (CURVEPOLYGON (LINESTRING (2 2, 3 2, 3 3, 2 3, 2 2)))");
    auto result = OverlayNG::overlay(a.get(), b.get(), OverlayNG::UNION);
    ensure(result != nullptr);
    ensure_not(result->getGeometryTypeId() == GeometryTypeId::GEOS_MULTIPOLYGON);
    ensure_not(result->getGeometryTypeId() == GeometryTypeId::GEOS_POLYGON);
    ensure(result->hasCurvedTypes());
}

} // namespace tut
