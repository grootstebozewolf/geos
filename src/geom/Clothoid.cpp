/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2026 Jeroen Bloemscheer
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/geom/Clothoid.h>

#include <geos/algorithm/Clothoid.h>
#include <geos/algorithm/CurveToLineParams.h>
#include <geos/constants.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/util.h>
#include <geos/util/IllegalArgumentException.h>

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace geos {
namespace geom {

namespace {

int
compareDouble(double a, double b)
{
    if (a < b) {
        return -1;
    }
    if (a > b) {
        return 1;
    }
    return 0;
}

std::unique_ptr<CoordinateSequence>
endpointsFor(const algorithm::clothoid::Parameters& parameters,
             bool hasZ, bool hasM, double z, double m)
{
    auto seq = std::make_unique<CoordinateSequence>(static_cast<std::size_t>(2), hasZ, hasM);
    double x1 = parameters.x0;
    double y1 = parameters.y0;
    if (parameters.length > 0.0) {
        algorithm::clothoid::position(parameters, parameters.length, x1, y1);
    }
    CoordinateXYZM start(parameters.x0, parameters.y0,
                         hasZ ? z : DoubleNotANumber,
                         hasM ? m : DoubleNotANumber);
    CoordinateXYZM end(x1, y1,
                       hasZ ? z : DoubleNotANumber,
                       hasM ? m : DoubleNotANumber);
    seq->setAt(start, 0);
    seq->setAt(end, 1);
    return seq;
}

} // namespace

Clothoid::Clothoid(const algorithm::clothoid::Parameters& parameters,
                   bool hasZ, bool hasM, double z, double m,
                   const GeometryFactory& factory)
    : SimpleCurve(endpointsFor(parameters, hasZ, hasM, z, m), false, factory),
      m_theta0(parameters.theta0),
      m_k0(parameters.k0),
      m_k1(parameters.k1),
      m_length(parameters.length)
{
    validateConstruction();
    envelope = computeClothoidEnvelope();
}

Clothoid::Clothoid(bool hasZ, bool hasM, const GeometryFactory& factory)
    : SimpleCurve(std::make_unique<CoordinateSequence>(static_cast<std::size_t>(0), hasZ, hasM), false, factory),
      m_theta0(0.0),
      m_k0(0.0),
      m_k1(0.0),
      m_length(0.0)
{
}

std::unique_ptr<Clothoid>
Clothoid::clone() const
{
    return std::unique_ptr<Clothoid>(cloneImpl());
}

std::string
Clothoid::getGeometryType() const
{
    return "Clothoid";
}

GeometryTypeId
Clothoid::getGeometryTypeId() const
{
    return GEOS_CLOTHOID;
}

double
Clothoid::getLength() const
{
    if (isEmpty()) {
        return 0.0;
    }
    return m_length;
}

double
Clothoid::getEndTangent() const
{
    return algorithm::clothoid::endHeading(parameters());
}

bool
Clothoid::hasCurvedComponents() const
{
    if (isEmpty() || !(m_length > 0.0)) {
        return false;
    }
    return m_k0 != 0.0 || m_k1 != 0.0;
}

void
Clothoid::normalize()
{
    if (isEmpty() || isClosed()) {
        return;
    }
    if (getEndCoordinate().compareTo(getStartCoordinate()) < 0) {
        std::unique_ptr<Clothoid> reversed(reverseImpl());
        points = reversed->points;
        envelope = reversed->envelope;
        m_theta0 = reversed->m_theta0;
        m_k0 = reversed->m_k0;
        m_k1 = reversed->m_k1;
        m_length = reversed->m_length;
    }
}

void
Clothoid::apply_rw(const CoordinateFilter*)
{
    throw util::UnsupportedOperationException(
        "Clothoid parameters are not edited through coordinate filters.");
}

void
Clothoid::apply_rw(CoordinateSequenceFilter&)
{
    throw util::UnsupportedOperationException(
        "Clothoid parameters are not edited through coordinate filters.");
}

bool
Clothoid::equalsExact(const Geometry* other, double tolerance) const
{
    if (!SimpleCurve::equalsExact(other, tolerance)) {
        return false;
    }
    const Clothoid* clothoid = detail::down_cast<const Clothoid*>(other);
    return std::abs(m_theta0 - clothoid->m_theta0) <= tolerance
        && std::abs(m_k0 - clothoid->m_k0) <= tolerance
        && std::abs(m_k1 - clothoid->m_k1) <= tolerance
        && std::abs(m_length - clothoid->m_length) <= tolerance;
}

bool
Clothoid::equalsIdentical(const Geometry* other) const
{
    if (!SimpleCurve::equalsIdentical(other)) {
        return false;
    }
    const auto& clothoid = static_cast<const Clothoid&>(*other);
    return m_theta0 == clothoid.m_theta0
        && m_k0 == clothoid.m_k0
        && m_k1 == clothoid.m_k1
        && m_length == clothoid.m_length;
}

int
Clothoid::compareToSameClass(const Geometry* geom) const
{
    const auto* clothoid = detail::down_cast<const Clothoid*>(geom);
    const int base = SimpleCurve::compareToSameClass(geom);
    if (base != 0) {
        return base;
    }
    int cmp = compareDouble(m_theta0, clothoid->m_theta0);
    if (cmp != 0) {
        return cmp;
    }
    cmp = compareDouble(m_k0, clothoid->m_k0);
    if (cmp != 0) {
        return cmp;
    }
    cmp = compareDouble(m_k1, clothoid->m_k1);
    if (cmp != 0) {
        return cmp;
    }
    return compareDouble(m_length, clothoid->m_length);
}

LineString*
Clothoid::getLinearizedImpl(const algorithm::CurveToLineParams& params) const
{
    if (isEmpty()) {
        return getFactory()->createLineString(hasZ(), hasM()).release();
    }

    const double stepRadians = params.getMaxStepDegrees() * MATH_PI / 180.0;
    std::size_t segments = algorithm::clothoid::linearizeSegmentCount(
        parameters(), stepRadians, params.getMaxDeviation());
    if (segments == 0) {
        segments = 1;
    }

    auto seq = std::make_unique<CoordinateSequence>(segments + 1, hasZ(), hasM());
    CoordinateXYZM ordinate;
    points->getAt(0, ordinate);
    const algorithm::clothoid::Parameters curve = parameters();
    for (std::size_t i = 0; i <= segments; ++i) {
        const double s = curve.length * (static_cast<double>(i) / static_cast<double>(segments));
        double x = curve.x0;
        double y = curve.y0;
        algorithm::clothoid::position(curve, s, x, y);
        ordinate.x = x;
        ordinate.y = y;
        seq->setAt(ordinate, i);
    }
    return getFactory()->createLineString(std::move(seq)).release();
}

Clothoid*
Clothoid::reverseImpl() const
{
    if (isEmpty()) {
        return cloneImpl();
    }

    algorithm::clothoid::Parameters reversed = parameters();
    double x1 = 0.0;
    double y1 = 0.0;
    algorithm::clothoid::position(reversed, reversed.length, x1, y1);
    reversed.x0 = x1;
    reversed.y0 = y1;
    reversed.theta0 = algorithm::clothoid::endHeading(parameters()) + MATH_PI;
    reversed.k0 = -m_k1;
    reversed.k1 = -m_k0;

    CoordinateXYZM start;
    points->getAt(0, start);
    return getFactory()->createClothoid(
        CoordinateXYZM(reversed.x0, reversed.y0, start.z, start.m),
        hasZ(), hasM(),
        reversed.theta0, reversed.k0, reversed.k1, reversed.length).release();
}

Envelope
Clothoid::computeClothoidEnvelope() const
{
    Envelope env;
    if (isEmpty()) {
        return env;
    }
    const algorithm::clothoid::Parameters curve = parameters();
    env.expandToInclude(curve.x0, curve.y0);
    double x1 = curve.x0;
    double y1 = curve.y0;
    algorithm::clothoid::position(curve, curve.length, x1, y1);
    env.expandToInclude(x1, y1);

    std::vector<double> stations;
    algorithm::clothoid::axisExtremaStations(curve, stations);
    for (double s : stations) {
        double x = 0.0;
        double y = 0.0;
        algorithm::clothoid::position(curve, s, x, y);
        env.expandToInclude(x, y);
    }
    return env;
}

void
Clothoid::validateConstruction() const
{
    if (!std::isfinite(m_theta0) || !std::isfinite(m_k0) || !std::isfinite(m_k1)
            || !std::isfinite(m_length) || !(m_length >= 0.0)) {
        throw util::IllegalArgumentException(
            "Clothoid tangent, curvatures, and length must be finite and length must be non-negative");
    }
    if (!std::isfinite(getStartCoordinate().x) || !std::isfinite(getStartCoordinate().y)) {
        throw util::IllegalArgumentException("Clothoid start coordinate must be finite");
    }
}

algorithm::clothoid::Parameters
Clothoid::parameters() const
{
    algorithm::clothoid::Parameters curve;
    curve.x0 = 0.0;
    curve.y0 = 0.0;
    curve.theta0 = m_theta0;
    curve.k0 = m_k0;
    curve.k1 = m_k1;
    curve.length = m_length;
    if (!isEmpty()) {
        curve.x0 = getStartCoordinate().x;
        curve.y0 = getStartCoordinate().y;
    }
    return curve;
}

}
}
