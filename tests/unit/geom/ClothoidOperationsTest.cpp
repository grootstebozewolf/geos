// Distance and point location with clothoids, alone and inside CompoundCurve,
// CurvePolygon and the Multi* collections. The reference is a linearization
// with a small maximum deviation of the clothoids.

#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>

#include <geos/algorithm/CurveToLineParams.h>
#include <geos/algorithm/PointLocator.h>
#include <geos/algorithm/locate/IndexedPointInAreaLocator.h>
#include <geos/algorithm/locate/SimplePointInAreaLocator.h>
#include <geos/constants.h>
#include <geos/geom/Clothoid.h>
#include <geos/geom/CompoundCurve.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CurvePolygon.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Location.h>
#include <geos/geom/MultiCurve.h>
#include <geos/geom/MultiSurface.h>
#include <geos/geom/Point.h>
#include <geos/io/WKTReader.h>
#include <geos/operation/distance/DistanceOp.h>
#include <geos/operation/distance/IndexedFacetDistance.h>

#include <cmath>
#include <string>
#include <vector>

using geos::algorithm::CurveToLineParams;
using geos::algorithm::PointLocator;
using geos::geom::Clothoid;
using geos::geom::ClothoidReferenceLocation;
using geos::geom::Coordinate;
using geos::geom::CoordinateXY;
using geos::geom::Geometry;
using geos::geom::Location;
using geos::geom::SimpleCurve;
using geos::operation::distance::DistanceOp;

namespace tut {

struct test_clothoidoperations_data {
    geos::geom::GeometryFactory::Ptr factory_ = geos::geom::GeometryFactory::create();
    geos::io::WKTReader reader_{*factory_};

    // Linearization deviation of the reference, and the agreement asked of
    // distances measured against it.
    static constexpr double DEVIATION = 1e-9;
    static constexpr double AGREEMENT = 1e-8;

    static ClothoidReferenceLocation at(double x, double y,
                                        double ux = 1.0, double uy = 0.0,
                                        double vx = 0.0, double vy = 1.0)
    {
        return ClothoidReferenceLocation{
            Coordinate(x, y, 0.0), Coordinate(ux, uy, 0.0), Coordinate(vx, vy, 0.0)};
    }

    std::unique_ptr<Clothoid> iso(const ClothoidReferenceLocation& loc,
                                  double scale, double t0, double t1) const
    {
        return factory_->createClothoid(loc, scale, t0, t1, false, false,
                                        geos::DoubleNotANumber, geos::DoubleNotANumber);
    }

    static std::unique_ptr<Geometry> linear(const Geometry& g)
    {
        return g.getLinearized(CurveToLineParams::maxDeviation(DEVIATION));
    }

    std::unique_ptr<Geometry> point(double x, double y) const
    {
        return factory_->createPoint(CoordinateXY(x, y));
    }

    /// The clothoid closed back to its start by a chord.
    std::unique_ptr<geos::geom::CompoundCurve> closed(std::unique_ptr<Clothoid> curve) const
    {
        const CoordinateXY start = curve->getStartCoordinate();
        const CoordinateXY end = curve->getEndCoordinate();
        std::vector<std::unique_ptr<SimpleCurve>> sections;
        sections.push_back(std::move(curve));
        sections.push_back(factory_->createLineString({end, start}));
        return factory_->createCompoundCurve(std::move(sections));
    }

    /// a holds the clothoids. b is measured exactly unless it holds a
    /// clothoid too, when it is also linearized.
    void ensureDistance(const std::string& what, const Geometry& a, const Geometry& b,
                        double deviation = DEVIATION, double agreement = AGREEMENT,
                        bool linearizeB = false) const
    {
        const auto params = CurveToLineParams::maxDeviation(deviation);
        const auto linearA = a.getLinearized(params);
        const auto linearB = linearizeB ? b.getLinearized(params) : b.clone();
        const double expected = linearA->distance(linearB.get());
        const double actual = DistanceOp::distance(a, b);
        ensure_equals(what, actual, expected, agreement);
        ensure_equals(what + " reversed", DistanceOp::distance(b, a), actual, agreement);
        // The nearest points are on the inputs and realize the distance.
        auto pts = DistanceOp::nearestPoints(&a, &b);
        const CoordinateXY p0 = pts->getAt<CoordinateXY>(0);
        const CoordinateXY p1 = pts->getAt<CoordinateXY>(1);
        ensure_equals(what + " nearest points", p0.distance(p1), actual, agreement);
        ensure(what + " first on a", linearA->distance(point(p0.x, p0.y).get()) <= agreement);
        ensure(what + " second on b", linearB->distance(point(p1.x, p1.y).get()) <= agreement);
    }
};

typedef test_group<test_clothoidoperations_data> group;
typedef group::object object;

group test_clothoidoperations_group("geos::geom::ClothoidOperations");

template<>
template<>
void object::test<1>()
{
    set_test_name("distance from points");

    auto spiral = iso(at(0.0, 0.0), 2.0, 0.0, 4.0);
    const double points[][2] = {
        {0.5, 1.0}, {1.0, 3.0}, {3.0, 0.0}, {-1.0, -1.0}, {2.5, 2.5}, {1.2, 1.4}, {10.0, 10.0}
    };
    for (const auto& xy : points) {
        ensureDistance("point " + std::to_string(xy[0]) + " " + std::to_string(xy[1]),
                       *spiral, *point(xy[0], xy[1]));
    }

    // Both sides of the inflection, sheared, and reflected.
    auto across = iso(at(1.0, -2.0, 2.0, 0.0, 0.5, 1.0), 1.5, -3.0, 2.5);
    auto reflected = iso(at(0.0, 0.0, 0.0, 1.0, 1.0, 0.0), 1.0, 2.0, -1.0);
    for (const auto& xy : points) {
        ensureDistance("across", *across, *point(xy[0], xy[1]));
        ensureDistance("reflected", *reflected, *point(xy[0], xy[1]));
    }
}

template<>
template<>
void object::test<2>()
{
    set_test_name("distance from segments, arcs and clothoids");

    auto spiral = iso(at(0.0, 0.0), 2.0, 0.0, 4.0);
    const char* others[] = {
        "LINESTRING (0 1, 3 1)",                              // crosses
        "LINESTRING (4 -1, 4 5)",                             // apart
        "LINESTRING (-1 -0.5, 3 -0.5, 3.5 3)",
        "CIRCULARSTRING (2 0, 3 1, 2 2)",
        "CIRCULARSTRING (-1 1, 0 2, 1 1)",
        "COMPOUNDCURVE ((5 0, 5 1), CIRCULARSTRING (5 1, 4.5 1.5, 5 2))",
    };
    for (const char* wkt : others) {
        ensureDistance(wkt, *spiral, *reader_.read(wkt));
    }
    ensure_equals("crossing", DistanceOp::distance(*spiral, *reader_.read("LINESTRING (0 1, 3 1)")), 0.0);

    // A clothoid pair is checked against a coarser reference: the dense
    // linear distance is quadratic in the vertex count.
    auto apart = iso(at(3.5, 3.0, 0.0, -1.0, 1.0, 0.0), 1.5, 0.5, 2.5);
    ensure(DistanceOp::distance(*spiral, *apart) > 0.1);
    ensureDistance("clothoid pair apart", *spiral, *apart, 1e-6, 1e-5, true);
    auto touching = iso(at(4.0, 0.0, -1.0, 0.0, 0.0, 1.0), 1.5, 0.5, 3.0);
    ensureDistance("clothoid pair touching", *spiral, *touching, 1e-6, 1e-5, true);
    auto crossing = iso(at(1.0, -1.0, 0.0, 1.0, -1.0, 0.0), 2.0, 0.0, 3.0);
    ensure(DistanceOp::distance(*spiral, *crossing) <= AGREEMENT);
}

template<>
template<>
void object::test<3>()
{
    set_test_name("point location on a clothoid");

    auto spiral = iso(at(0.0, 0.0), 2.0, 0.0, 4.0);
    auto line = linear(*spiral);
    PointLocator locator;
    // Linearization vertices are evaluated points of the curve.
    const auto seq = line->getCoordinates();
    for (std::size_t i = 1; i + 1 < seq->size(); i += 97) {
        ensure_equals(locator.locate(seq->getAt<CoordinateXY>(i), spiral.get()), Location::INTERIOR);
    }
    ensure_equals(locator.locate(spiral->getStartCoordinate(), spiral.get()), Location::BOUNDARY);
    ensure_equals(locator.locate(spiral->getEndCoordinate(), spiral.get()), Location::BOUNDARY);
    ensure_equals(locator.locate(CoordinateXY(1.0, 1.0), spiral.get()), Location::EXTERIOR);

    std::vector<std::unique_ptr<SimpleCurve>> sections;
    const CoordinateXY end = spiral->getEndCoordinate();
    sections.push_back(spiral->clone());
    sections.push_back(factory_->createLineString({end, CoordinateXY(end.x + 1.0, end.y)}));
    auto compound = factory_->createCompoundCurve(std::move(sections));
    const CoordinateXY onClothoid = seq->getAt<CoordinateXY>(50);
    ensure_equals(locator.locate(onClothoid, static_cast<const Geometry*>(compound.get())), Location::INTERIOR);
}

template<>
template<>
void object::test<4>()
{
    set_test_name("point in a CurvePolygon with a clothoid ring section");

    using geos::algorithm::locate::SimplePointInAreaLocator;
    auto polygon = factory_->createCurvePolygon(closed(iso(at(0.0, 0.0), 2.0, 0.0, 3.0)));
    auto reference = linear(*polygon);
    PointLocator locator;
    const auto& env = *polygon->getEnvelopeInternal();
    int inside = 0;
    for (int i = 0; i <= 40; i++) {
        for (int j = 0; j <= 40; j++) {
            const CoordinateXY p(env.getMinX() + env.getWidth() * i / 40.0,
                                 env.getMinY() + env.getHeight() * j / 40.0);
            auto pt = point(p.x, p.y);
            if (reference->getBoundary()->distance(pt.get()) < 1e-6) {
                continue;
            }
            const Location expected = SimplePointInAreaLocator::locate(p, reference.get());
            ensure_equals(SimplePointInAreaLocator::locate(p, polygon.get()), expected);
            ensure_equals(locator.locate(p, static_cast<const Geometry*>(polygon.get())), expected);
            inside += (expected == Location::INTERIOR);
        }
    }
    ensure(inside > 50);

    // On the clothoid section.
    const CoordinateXY onRing = linear(*polygon->getExteriorRing())->getCoordinates()->getAt<CoordinateXY>(100);
    ensure_equals(SimplePointInAreaLocator::locate(onRing, polygon.get()), Location::BOUNDARY);
}

template<>
template<>
void object::test<5>()
{
    set_test_name("Multi* and collections");

    auto spiral = iso(at(0.0, 0.0), 2.0, 0.0, 3.0);
    std::vector<std::unique_ptr<Geometry>> curves;
    curves.push_back(spiral->clone());
    curves.push_back(reader_.read("CIRCULARSTRING (4 0, 5 1, 6 0)"));
    auto multiCurve = factory_->createMultiCurve(std::move(curves));
    ensureDistance("MultiCurve to point", *multiCurve, *point(2.0, 2.0));
    ensureDistance("MultiCurve to line", *multiCurve, *reader_.read("LINESTRING (3 3, 3.5 2)"));

    std::vector<std::unique_ptr<Geometry>> surfaces;
    surfaces.push_back(factory_->createCurvePolygon(closed(spiral->clone())));
    surfaces.push_back(reader_.read("POLYGON ((10 0, 11 0, 11 1, 10 1, 10 0))"));
    auto multiSurface = factory_->createMultiSurface(std::move(surfaces));
    PointLocator locator;
    ensure_equals(locator.locate(CoordinateXY(1.5, 0.4), multiSurface.get()), Location::INTERIOR);
    ensure_equals(locator.locate(CoordinateXY(1.0, 1.5), multiSurface.get()), Location::EXTERIOR);
    ensure_equals(locator.locate(CoordinateXY(10.5, 0.5), multiSurface.get()), Location::INTERIOR);

    // A point inside the surface is at distance zero from it.
    ensure_equals(DistanceOp::distance(*multiSurface, *point(1.5, 0.4)), 0.0);
    ensureDistance("MultiSurface to point", *multiSurface, *point(1.0, 2.0));
}

template<>
template<>
void object::test<6>()
{
    set_test_name("large turn across the inflection");

    // Each side turns 16200 rad. A fine reference would need tens of
    // millions of points, so this one is coarser.
    auto spiral = iso(at(0.0, 0.0, 2.0, 0.0, 0.5, 1.0), 1.0, -180.0, 180.0);
    auto reference = spiral->getLinearized(CurveToLineParams::maxDeviation(1e-5));
    auto pt = point(0.3, 0.2);
    ensure_equals("point", DistanceOp::distance(*spiral, *pt), reference->distance(pt.get()), 1e-4);
    auto line = reader_.read("LINESTRING (-3 3, 3 3)");
    ensure_equals("segment", DistanceOp::distance(*spiral, *line), reference->distance(line.get()), 1e-4);
}

template<>
template<>
void object::test<7>()
{
    set_test_name("clothoids crossing away from their vertices touch at distance zero");

    // An S-shaped spiral through (0, 0.5), and its image rotated a quarter
    // turn about (0.3, 0.5). They cross at a point that is not a vertex of
    // either, which the pair search finds only to within rounding.
    auto a = iso(at(0.0, 0.5), 2.0, -2.0, 2.0);
    auto b = iso(at(0.3, 0.5, 0.0, 1.0, -1.0, 0.0), 2.0, -2.0, 2.0);

    ensure_equals("distance", DistanceOp::distance(*a, *b), 0.0);
    ensure_equals("distance reversed", DistanceOp::distance(*b, *a), 0.0);
    ensure("intersects", a->intersects(b.get()));
    ensure("intersects reversed", b->intersects(a.get()));
    ensure("linearized inputs intersect", linear(*a)->intersects(linear(*b).get()));
}

template<>
template<>
void object::test<8>()
{
    set_test_name("indexed point location in a ring with a clothoid section");

    using geos::algorithm::locate::IndexedPointInAreaLocator;
    using geos::algorithm::locate::SimplePointInAreaLocator;

    auto spiral = iso(at(0.0, 0.5), 2.0, -2.0, 2.0);
    const CoordinateXY start = spiral->getStartCoordinate();
    const CoordinateXY end = spiral->getEndCoordinate();
    std::vector<std::unique_ptr<SimpleCurve>> sections;
    sections.push_back(std::move(spiral));
    sections.push_back(factory_->createLineString({end, CoordinateXY(2, -2), CoordinateXY(-2, -2), start}));
    auto polygon = factory_->createCurvePolygon(factory_->createCompoundCurve(std::move(sections)));

    IndexedPointInAreaLocator locator(*polygon);
    auto check = [&](double x, double y, Location expected) {
        const CoordinateXY q(x, y);
        ensure_equals("indexed", locator.locate(&q), expected);
        ensure_equals("simple", SimplePointInAreaLocator::locate(q, polygon.get()), expected);
    };
    check(0.0, 0.5, Location::BOUNDARY);  // inflection point of the spiral
    check(0.0, 0.4, Location::INTERIOR);
    check(0.0, -2.0, Location::BOUNDARY);
    // The chord runs through (-1, 0.33) and (1, 0.67), and the spiral through
    // about (-1, 0.46) and (1, 0.54), so these sit between the two.
    check(1.0, 0.6, Location::EXTERIOR);
    check(-1.0, 0.4, Location::INTERIOR);
}

template<>
template<>
void object::test<9>()
{
    set_test_name("indexed facet distance to clothoids");

    using geos::operation::distance::IndexedFacetDistance;

    auto spiral = iso(at(0.0, 0.5), 2.0, -2.0, 2.0);

    // Directly below the inflection, where the spiral is level: 0.1, where
    // its chord would give less.
    ensure_equals("point", IndexedFacetDistance::distance(spiral.get(), point(0.0, 0.4).get()), 0.1, 1e-12);

    std::vector<std::string> targets = {
        "LINESTRING (-3 2, 3 2.5)",
        "CIRCULARSTRING (-1 -1, 0 -0.5, 1 -1)",
        "LINESTRING (0 -1, 0 2)",
        "MULTICURVE (CIRCULARSTRING (3 0, 4 1, 5 0), (6 6, 7 7))",
    };
    for (const auto& wkt : targets) {
        auto target = reader_.read(wkt);
        const double expected = DistanceOp::distance(*spiral, *target);
        ensure_equals(wkt, IndexedFacetDistance::distance(spiral.get(), target.get()), expected, 1e-12);
        ensure_equals(wkt + " reversed", IndexedFacetDistance::distance(target.get(), spiral.get()), expected, 1e-12);
        auto pts = IndexedFacetDistance::nearestPoints(spiral.get(), target.get());
        const double nearest = pts->getAt<CoordinateXY>(0).distance(pts->getAt<CoordinateXY>(1));
        ensure_equals(wkt + " nearest points", nearest, expected, 1e-12);
    }

    auto other = iso(at(0.3, 0.5, 0.0, 1.0, -1.0, 0.0), 2.0, -2.0, 2.0);
    ensure_equals("crossing clothoid", IndexedFacetDistance::distance(spiral.get(), other.get()), 0.0);
}

} // namespace tut
