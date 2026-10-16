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

#pragma once

#include <geos/geom/LineString.h>
#include <geos/geom/SimpleCurve.h>

namespace geos {
namespace algorithm {
namespace clothoid {
struct Parameters;
}
}
namespace geom {

/**
 * A single clothoid (Euler spiral) segment.
 *
 * Provisional fork extension. The WKT keyword CLOTHOID and WKB type code 18
 * are not an OGC or ISO/IEC 13249-3 assignment. This type stores one
 * arc-length segment
 *
 *   CLOTHOID [Z][M] (x y [z] [m], theta0, k0, k1, length)
 *
 * with theta0 in radians from +X, k0 and k1 the endpoint curvatures, and
 * length the arc length. Z and M, when present, are constant along the
 * segment in this slice.
 *
 * The coordinate sequence holds the start point and the integrated end
 * point. It is not a sampling of the spiral. getLinearized builds that
 * sampling. Topology operations that refuse curved components still refuse
 * a clothoid whose curvature is not identically zero; this type does not
 * implement a clothoid DE-9IM, overlay, or buffer.
 *
 * allure-wkt's COMPOUNDCURVE member form CLOTHOID (k0, k1, length), which
 * inherits its start point and tangent, is not accepted here.
 */
class GEOS_DLL Clothoid : public SimpleCurve {

public:

    friend class GeometryFactory;

    ~Clothoid() override = default;

    std::unique_ptr<Clothoid> clone() const;

    std::string getGeometryType() const override;

    GeometryTypeId getGeometryTypeId() const override;

    double getLength() const override;

    double getStartTangent() const
    {
        return m_theta0;
    }

    double getStartCurvature() const
    {
        return m_k0;
    }

    double getEndCurvature() const
    {
        return m_k1;
    }

    /// theta(length). Equals getStartTangent when the clothoid is empty.
    double getEndTangent() const;

    bool hasCurvedComponents() const override;

    bool hasCurvedTypes() const override
    {
        return true;
    }

    bool isCurved() const override
    {
        return hasCurvedComponents();
    }

    void normalize() override;

    std::unique_ptr<Clothoid> reverse() const
    {
        return std::unique_ptr<Clothoid>(reverseImpl());
    }

    void apply_rw(const CoordinateFilter* filter) override;

    void apply_rw(CoordinateSequenceFilter& filter) override;

    bool equalsExact(const Geometry* other, double tolerance = 0) const override;

    bool equalsIdentical(const Geometry* other) const override;

protected:

    Clothoid(const algorithm::clothoid::Parameters& parameters,
             bool hasZ, bool hasM, double z, double m,
             const GeometryFactory& factory);

    Clothoid(bool hasZ, bool hasM, const GeometryFactory& factory);

    Clothoid(const Clothoid& other) = default;

    Clothoid* cloneImpl() const override
    {
        return new Clothoid(*this);
    }

    int compareToSameClass(const Geometry* geom) const override;

    void geometryChangedAction() override
    {
        envelope = computeClothoidEnvelope();
    }

    Clothoid* getCurvedImpl(const algorithm::LineToCurveParams&) const override
    {
        return cloneImpl();
    }

    LineString* getLinearizedImpl(const algorithm::CurveToLineParams& params) const override;

    int getSortIndex() const override
    {
        return SORTINDEX_CLOTHOID;
    }

    Clothoid* reverseImpl() const override;

private:

    Envelope computeClothoidEnvelope() const;

    void validateConstruction() const;

    algorithm::clothoid::Parameters parameters() const;

    double m_theta0;
    double m_k0;
    double m_k1;
    double m_length;
};

}
}
