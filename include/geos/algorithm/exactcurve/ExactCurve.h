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

#pragma once

#include <geos/export.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Geometry.h>

#include <memory>

namespace geos {
namespace algorithm {
namespace exactcurve {

/**
 * Thin exact-curve protocol: six operations, no densify on
 * length / pointAt / isExact. toLinear is the only densify path.
 */
class GEOS_DLL ExactCurve {
public:
    virtual ~ExactCurve() = default;
    virtual const geom::CoordinateXY& getStart() const = 0;
    virtual const geom::CoordinateXY& getEnd() const = 0;
    virtual double length() const = 0;
    virtual geom::CoordinateXY pointAt(double t) const = 0; // t ∈ [0,1]
    virtual std::unique_ptr<geom::Geometry> toLinear(double tolerance) const = 0;
    virtual bool isExact() const = 0;
};

} // namespace exactcurve
} // namespace algorithm
} // namespace geos
