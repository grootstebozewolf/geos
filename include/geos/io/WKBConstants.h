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
        // Provisional fork extension. Not an OGC, ISO, or SQL/MM type code.
        // Body layout (after the standard byte order, type word, and optional SRID):
        //   x, y, [z], [m], theta0, k0, k1, length
        // Empty clothoids are encoded with NaN ordinates, matching POINT EMPTY.
        wkbClothoid = 18,
    };

    enum wkbFlavour {
        wkbExtended = 1,
        wkbIso = 2
    };

}

} // namespace geos::io
} // namespace geos

