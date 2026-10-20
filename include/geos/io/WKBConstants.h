/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2005-2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: io/WKBConstants.java rev. 1.1 (JTS-1.10)
 *
 **********************************************************************/

#pragma once

namespace geos {
namespace io {

/// Constant values used by the WKB format
namespace WKBConstants {

    enum byteOrder {
        wkbXDR = 0,
        wkbNDR = 1
    };

    enum wkbType {
        wkbPoint = 1,
        wkbLineString = 2,
        wkbPolygon = 3,
        wkbMultiPoint = 4,
        wkbMultiLineString = 5,
        wkbMultiPolygon = 6,
        wkbGeometryCollection = 7,
        wkbCircularString = 8,
        wkbCompoundCurve = 9,
        wkbCurvePolygon = 10,
        wkbMultiCurve = 11,
        wkbMultiSurface = 12,
        // ISO/IEC 13249-3:2016 §5.1.68 Table 15. ST_Circle is not implemented.
        // 18 is not a clothoid; ST_Clothoid is 22.
        wkbCircle = 18,
        // ISO/IEC 13249-3:2016 §5.1.68 Table 15 <wkbclothoid>.
        // Z / M / ZM are 1022 / 2022 / 3022. The body is an affine placement
        // (§5.1.68 <wkbreferencelocation>), scale factor, and start/end distance.
        // Empty is the type word with that body omitted (§5.1.68 bg).
        wkbClothoid = 22,
        // ISO/IEC 13249-3:2016 §5.1.68 Table 15. ST_SpiralCurve is not implemented.
        wkbSpiralCurve = 23,
        // Support type embedded in clothoid WKB. Z form is 1102.
        wkbAffinePlacement = 102,
    };

    enum wkbFlavour {
        wkbExtended = 1,
        wkbIso = 2
    };

}

} // namespace geos::io
} // namespace geos

