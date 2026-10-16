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

#include <geos/export.h>

#include <cstddef>
#include <vector>

namespace geos {
namespace algorithm {
namespace clothoid {

/**
 * Arc-length clothoid with linear curvature.
 *
 * Heading and curvature on s in [0, length], when length > 0:
 *   kappa(s) = k0 + (k1 - k0) * (s / length)
 *   theta(s) = theta0 + k0 * s + 0.5 * (k1 - k0) / length * s^2
 *
 * Position is the unit-speed integral of (cos theta, sin theta).
 * When the curvature rate is nonzero and the completed-square shift stays
 * moderate, that integral is reduced to a Fresnel increment and evaluated
 * numerically. Constant curvature uses the circular (or straight)
 * antiderivative. Other cases use parameter quadrature.
 *
 * These routines evaluate the integrals. They do not decide a DE-9IM
 * matrix, an overlay, or a buffer.
 */
struct GEOS_DLL Parameters {
    double x0;
    double y0;
    double theta0;
    double k0;
    double k1;
    double length;
};

struct GEOS_DLL Moments {
    double P;
    double Q;
    double R;
    double T;
    double S2c;
    double S2s;
};

struct GEOS_DLL LengthSolve {
    double length;
    int iterations;
    bool converged;
};

/// Default residual scale used by the reference Halley loop.
constexpr double LENGTH_TOLERANCE = 1e-13;
constexpr int LENGTH_MAX_ITERATIONS = 50;

/**
 * Fresnel integrals of ClothoidFresnel.v:
 *   fresnelC2(x) = ∫_0^x cos(t^2 / 2) dt
 *   fresnelS2(x) = ∫_0^x sin(t^2 / 2) dt
 *
 * Power series on [-5, 5]. Outside that window the increment is a
 * 16-point Gauss–Legendre panel sum of width 0.5. Both are odd.
 */
GEOS_DLL double fresnelC2(double x);
GEOS_DLL double fresnelS2(double x);

/// Heading theta(s). Outside [0, length] this is the same polynomial.
GEOS_DLL double heading(const Parameters& curve, double s);

/// Curvature kappa(s).
GEOS_DLL double curvature(const Parameters& curve, double s);

/// theta(length) = theta0 + 0.5 * (k0 + k1) * length, when length > 0.
GEOS_DLL double endHeading(const Parameters& curve);

/// Unit-speed position at station s.
GEOS_DLL void position(const Parameters& curve, double s, double& x, double& y);

/**
 * Six moment integrals on tau in [0, 1],
 *   psi(tau) = k0 * tau + 0.5 * (k1 - k0) * tau^2
 *   (P, Q, R, T, S2c, S2s) = ∫ {1, psi, psi^2} * {cos, sin}(L * psi) d tau
 * using the 32-point Gauss–Legendre rule from the Halley reference solver.
 */
GEOS_DLL Moments moments(double length, double k0, double k1);

/**
 * Halley solve of f(L) = L^2 (P^2 + Q^2) - d^2 for the chord from
 * (x0, y0) to (x1, y1). converged is true when the reference residual
 * test passes inside the iteration budget. A true flag is not a claim
 * that the root is the only positive root.
 */
GEOS_DLL LengthSolve solveLengthHalley(
    double x0, double y0, double x1, double y1,
    double k0, double k1,
    double tolerance = LENGTH_TOLERANCE,
    int maxIterations = LENGTH_MAX_ITERATIONS);

/// Same residual, Newton update. Included so the reference pair stays together.
GEOS_DLL LengthSolve solveLengthNewton(
    double x0, double y0, double x1, double y1,
    double k0, double k1,
    double tolerance = LENGTH_TOLERANCE,
    int maxIterations = LENGTH_MAX_ITERATIONS);

/**
 * Place a clothoid on the chord by solving L, then rotating the local
 * frame (theta0 = 0) onto the chord direction. When converged is false
 * the parameters are the last iterate and are not claimed to meet the end point.
 */
struct GEOS_DLL ChordPlacement {
    Parameters parameters;
    LengthSolve solve;
};

GEOS_DLL ChordPlacement placeOnChord(
    double x0, double y0, double x1, double y1,
    double k0, double k1,
    double tolerance = LENGTH_TOLERANCE,
    int maxIterations = LENGTH_MAX_ITERATIONS);

/**
 * Stations in (0, length) where theta(s) is an integer multiple of pi/2.
 * Those are the interior points ClothoidEnvelope.v uses for the axis-aligned box.
 * Endpoints are not included.
 */
GEOS_DLL void axisExtremaStations(const Parameters& curve, std::vector<double>& stations);

/**
 * Segment count for a uniform arc-length sample.
 *
 * ClothoidLinearize.v uses n >= 2 and tol = kappa_max * h^2 / 8,
 * h = length / n, kappa_max = max(|k0|, |k1|).
 * maxStepRadians limits the heading change per segment by kappa_max * h.
 * A straight curve (kappa_max == 0) returns 1. length <= 0 returns 0.
 */
GEOS_DLL std::size_t linearizeSegmentCount(
    const Parameters& curve,
    double maxStepRadians,
    double maxDeviation);

}
}
}
