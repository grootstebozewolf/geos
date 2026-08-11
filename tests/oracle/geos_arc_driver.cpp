// GEOS-side differential driver speaking the RocqRefRunner (oracle_bin)
// stdin/stdout protocol for the arc modes, so the NetTopologySuite.Proofs
// harnesses can gate GEOS's CircularArcIntersector against the verified
// oracle and its exact-rational ground truth.
//
// Supported modes:
//   ARC_ARC_XY      6 lines "x y" (arc1 A B C, arc2 A B C)
//                   -> "N [x y]*" | "DEGENERATE" | "COINCIDENT" | "NAN"
//   ARC_SEGMENT_XY  5 lines "x y" (arc A B C, segment P Q)
//                   -> "N [x y]*" | "DEGENERATE" | "NAN"
//
// Coordinates are printed as C hex floats (%a) to match the oracle's %h.
// Modes loop until EOF, one result line per case.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include <geos/algorithm/CircularArcIntersector.h>
#include <geos/geom/CircularArc.h>
#include <geos/geom/CoordinateSequence.h>

using geos::algorithm::CircularArcIntersector;
using geos::geom::CircularArc;
using geos::geom::CoordinateSequence;
using geos::geom::CoordinateXY;

static bool
readPoint(CoordinateXY& pt)
{
    std::string line;
    if (!std::getline(std::cin, line)) {
        return false;
    }
    // strtod accepts both decimal and C99 hex-float ("0x1.8p+1") forms.
    const char* s = line.c_str();
    char* end = nullptr;
    pt.x = std::strtod(s, &end);
    pt.y = std::strtod(end, &end);
    return true;
}

static bool
finite(const CoordinateXY& pt)
{
    return std::isfinite(pt.x) && std::isfinite(pt.y);
}

static void
printPoints(const CircularArcIntersector& cai)
{
    std::printf("%d", cai.getNumPoints());
    for (std::uint8_t i = 0; i < cai.getNumPoints(); i++) {
        const auto& pt = cai.getPoint(i);
        std::printf(" %a %a", pt.x, pt.y);
    }
    std::printf("\n");
}

static void
runArcArc()
{
    CoordinateXY a1, b1, c1, a2, b2, c2;
    if (!(readPoint(a1) && readPoint(b1) && readPoint(c1) &&
          readPoint(a2) && readPoint(b2) && readPoint(c2))) {
        std::exit(1);
    }

    if (!(finite(a1) && finite(b1) && finite(c1) &&
          finite(a2) && finite(b2) && finite(c2))) {
        std::printf("NAN\n");
        return;
    }

    CircularArc arc1 = CircularArc::create(a1, b1, c1);
    CircularArc arc2 = CircularArc::create(a2, b2, c2);

    // The oracle reports DEGENERATE when either arc's controls are collinear.
    if (arc1.isLinear() || arc2.isLinear()) {
        std::printf("DEGENERATE\n");
        return;
    }

    CircularArcIntersector cai;
    cai.intersects(arc1, arc2);

    if (cai.getResult() == CircularArcIntersector::COCIRCULAR_INTERSECTION ||
        cai.getNumArcs() > 0) {
        std::printf("COINCIDENT nPt=%d nArc=%d\n", cai.getNumPoints(), cai.getNumArcs());
        return;
    }

    printPoints(cai);
}

static void
runArcSegment()
{
    CoordinateXY a, b, c, p, q;
    if (!(readPoint(a) && readPoint(b) && readPoint(c) &&
          readPoint(p) && readPoint(q))) {
        std::exit(1);
    }

    if (!(finite(a) && finite(b) && finite(c) && finite(p) && finite(q))) {
        std::printf("NAN\n");
        return;
    }

    CircularArc arc = CircularArc::create(a, b, c);

    // Match the oracle: collinear arc OR zero-length segment is DEGENERATE.
    if (arc.isLinear() || (p.x == q.x && p.y == q.y)) {
        std::printf("DEGENERATE\n");
        return;
    }

    CoordinateSequence seq(2, false, false);
    seq.setAt(p, 0);
    seq.setAt(q, 1);

    CircularArcIntersector cai;
    cai.intersects(arc, seq, 0, 1, false);

    printPoints(cai);
}

int
main()
{
    std::string mode;
    while (std::getline(std::cin, mode)) {
        // Skip blank separator lines between stanzas.
        if (mode.find_first_not_of(" \t\r") == std::string::npos) {
            continue;
        }
        if (!mode.empty() && mode.back() == '\r') {
            mode.pop_back();
        }

        if (mode == "ARC_ARC_XY") {
            runArcArc();
        } else if (mode == "ARC_SEGMENT_XY") {
            runArcSegment();
        } else {
            std::fprintf(stderr, "geos_arc_driver: unknown mode: %s\n", mode.c_str());
            return 1;
        }
        std::fflush(stdout);
    }
    return 0;
}
