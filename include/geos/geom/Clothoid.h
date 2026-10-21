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
 * Reference location of an ST_Clothoid.
 *
 * ISO/IEC 13249-3:2016 §4.9.3. Two reference directions, the images of the
 * constructive X and Y axes. §7.8.1 requires that cardinality to be 2.
 * Z components are zero when the clothoid is not 3D.
 */
struct GEOS_DLL ClothoidReferenceLocation {
    double x;
    double y;
    double z;
    double ux;
    double uy;
    double uz;
    double vx;
    double vy;
    double vz;
};

/**
 * A single ST_Clothoid segment.
 *
 * ISO/IEC 13249-3:2016 §4.2.11 and §7.8. WKT keyword CLOTHOID is the SQL/MM
 * keyword (§5.1.67). WKB type 22, with Z/M/ZM codes 1022/2022/3022
 * (§5.1.68 Table 15). Type 18 is ST_Circle in that table and is not a clothoid.
 *
 * The stored attributes are the reference location (an ST_AffinePlacement),
 * the scale factor, and the start and end constructive distances. Optional
 * measure is a start/end pair (§4.2.11, §7.8.1), not an ordinate of the
 * location point.
 *
 * In the constructive frame the curve is the unit-speed Fresnel clothoid
 * whose scale factor A satisfies A*A = R*t (ISO 19107 clothoid, the curve
 * SQL/MM §4.2.11 places with the Fresnel integrals):
 *
 *   theta(t) = t^2 / (2 A^2)
 *   kappa(t) = t / A^2
 *   X(t) = integral_0^t cos(tau^2 / (2 A^2)) d tau
 *   Y(t) = integral_0^t sin(tau^2 / (2 A^2)) d tau
 *
 * t is the arc length from the inflection in that frame, and the Fresnel
 * limit. The affine placement maps (X, Y) into the coordinate system:
 *
 *   P = location + X * refDirection[0] + Y * refDirection[1]
 *
 * A negative scale factor is stored so WKT and WKB round-trip. ISO/IEC
 * 13249-3:2016 §7.8.4 observes and mutates ST_ScaleFactor and states no sign
 * rule. Evaluation uses A*A, so a negative scale currently matches |A|.
 * WKT and WKB accept only this ISO body. GeometryFactory::createClothoid
 * still accepts an arc-length pose (start, theta0, k0, k1, length) and
 * converts it; that tuple is not a CLOTHOID encoding.
 *
 * The coordinate sequence holds the start point and the integrated end
 * point. getLinearized builds the sampling. Topology operations that refuse
 * curved components still refuse a clothoid whose curvature is not
 * identically zero. This type does not implement a clothoid DE-9IM, overlay,
 * or buffer. ST_SpiralCurve (§4.2.12, §7.9, WKB 23) is not implemented.
 * allure-wkt's inherited-start form CLOTHOID (k0, k1, length) is not accepted.
 */
class GEOS_DLL Clothoid : public SimpleCurve {

public:

    friend class GeometryFactory;

    ~Clothoid() override = default;

    std::unique_ptr<Clothoid> clone() const;

    std::string getGeometryType() const override;

    GeometryTypeId getGeometryTypeId() const override;

    double getLength() const override;

    /// Travel heading at the start, radians from +X. Zero when empty.
    double getStartTangent() const;

    /// Signed curvature at the start. Zero when empty.
    double getStartCurvature() const;

    /// Signed curvature at the end. Zero when empty.
    double getEndCurvature() const;

    /// theta at the end of travel. Equals getStartTangent when empty.
    double getEndTangent() const;

    /// ISO/IEC 13249-3:2016 §7.8 ST_ScaleFactor. Zero when empty.
    double getScaleFactor() const
    {
        return m_scale;
    }

    /// Constructive distance of the start from the inflection (§7.8).
    double getStartDistance() const
    {
        return m_t0;
    }

    /// Constructive distance of the end from the inflection (§7.8).
    double getEndDistance() const
    {
        return m_t1;
    }

    /// ST_StartM. Zero when the clothoid is not measured.
    double getStartM() const
    {
        return m_startM;
    }

    /// ST_EndM. Zero when the clothoid is not measured.
    double getEndM() const
    {
        return m_endM;
    }

    ClothoidReferenceLocation getReferenceLocation() const;

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

    Clothoid(const ClothoidReferenceLocation& location,
             double scaleFactor,
             double startDistance,
             double endDistance,
             bool hasZ, bool hasM,
             double startM, double endM,
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

    void loadConstructive(algorithm::clothoid::Parameters& curve) const;

    void pointOnCurve(double travel, double& x, double& y, double& z, double& m) const;

    Envelope computeClothoidEnvelope() const;

    void validateConstruction() const;

    double travelHeading(double constructiveDistance) const;

    double travelCurvature(double constructiveDistance) const;

    double m_locX;
    double m_locY;
    double m_locZ;
    double m_ux;
    double m_uy;
    double m_uz;
    double m_vx;
    double m_vy;
    double m_vz;
    double m_scale;
    double m_t0;
    double m_t1;
    double m_startM;
    double m_endM;
};

}
}
