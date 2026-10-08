// geos_probe: batch driver that exposes GEOS circular-curve primitives to the
// differential bug hunter (hunt.py).  One case per stdin line, one JSON object
// per stdout line.  Numbers are printed with 17 significant digits so that the
// Python side sees the exact binary64 values GEOS produced.
//
// Case grammar (whitespace separated; coordinates are decimal or %a hex floats):
//   A  x0 y0 x1 y1 x2 y2 dev            single arc properties + linearization
//   AA a0x a0y a1x a1y a2x a2y b0x b0y b1x b1y b2x b2y   arc/arc
//   AS a0x a0y a1x a1y a2x a2y s0x s0y s1x s1y           arc/segment
//   AP a0x a0y a1x a1y a2x a2y px py                     arc/point
//   SS p0x p0y p1x p1y q0x q0y q1x q1y                   segment/segment via CircularArcIntersector
//   G  <WKT>                                             length, area, linearize
//   N  <0|1> <WKT>                                       GeometryNoder (flag = preserveCompoundCurves)
#include <geos/algorithm/CircularArcIntersector.h>
#include <geos/algorithm/CurveToLineParams.h>
#include <geos/geom/CircularArc.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Curve.h>
#include <geos/geom/SimpleCurve.h>
#include <geos/io/WKTReader.h>
#include <geos/noding/GeometryNoder.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace geos;
using geom::CircularArc;
using geom::CoordinateSequence;
using geom::CoordinateXY;
using algorithm::CircularArcIntersector;

static std::string num(double d) {
    if (std::isnan(d)) return "\"nan\"";
    if (std::isinf(d)) return d > 0 ? "\"inf\"" : "\"-inf\"";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", d);
    return buf;
}
static std::string pt(const CoordinateXY& p) { return "[" + num(p.x) + "," + num(p.y) + "]"; }

static double rd(std::istream& in) {
    std::string s; in >> s;
    if (s.empty()) throw std::runtime_error("missing number");
    return std::strtod(s.c_str(), nullptr);
}
static CoordinateXY rp(std::istream& in) { double x = rd(in); double y = rd(in); return {x, y}; }

static std::unique_ptr<CoordinateSequence> seq3(const CoordinateXY& a, const CoordinateXY& b, const CoordinateXY& c) {
    auto s = std::make_unique<CoordinateSequence>(3u, false, false);
    s->setAt(a, 0); s->setAt(b, 1); s->setAt(c, 2);
    return s;
}

static std::string arcProps(const CircularArc& arc, double dev) {
    std::ostringstream o;
    o << "\"orient\":" << arc.getOrientation()
      << ",\"center\":" << pt(arc.getCenter())
      << ",\"radius\":" << num(arc.getRadius())
      << ",\"angle\":" << num(arc.getAngle())
      << ",\"length\":" << num(arc.getLength())
      << ",\"area\":" << num(arc.getArea())
      << ",\"isCircle\":" << (arc.isCircle() ? "true" : "false")
      << ",\"isLinear\":" << (arc.isLinear() ? "true" : "false");
    if (dev > 0) {
        CoordinateSequence lin(0u, false, false);
        auto params = algorithm::CurveToLineParams::maxDeviation(dev);
        arc.addLinearizedPoints(lin, params);
        o << ",\"step\":" << num(params.getStepSizeDegrees(arc)) << ",\"nlin\":" << lin.size() << ",\"lin\":[";
        const std::size_t cap = 100000;
        for (std::size_t i = 0; i < lin.size(); i++) {
            if (lin.size() > cap && i > 10 && i + 10 < lin.size()) continue;   // huge output: keep ends only
            if (i) o << ",";
            o << pt(lin.getAt<CoordinateXY>(i));
        }
        o << "]";
        if (lin.size() > cap) o << ",\"linTruncated\":true";
    }
    return o.str();
}

static std::string ciResult(const CircularArcIntersector& ci) {
    std::ostringstream o;
    o << "{\"type\":" << int(ci.getResult()) << ",\"nPt\":" << int(ci.getNumPoints())
      << ",\"nArc\":" << int(ci.getNumArcs()) << ",\"pts\":[";
    for (int i = 0; i < ci.getNumPoints() && i < 2; i++) {
        if (i) o << ",";
        o << pt(ci.getPoint(static_cast<std::uint8_t>(i)));
    }
    o << "],\"arcs\":[";
    for (int i = 0; i < ci.getNumArcs() && i < 2; i++) {
        if (i) o << ",";
        const auto& a = ci.getArc(static_cast<std::uint8_t>(i));
        o << "[" << pt(a.p0()) << "," << pt(a.p1()) << "," << pt(a.p2()) << "]";
    }
    o << "]}";
    return o.str();
}

static void dumpGeom(const geom::Geometry& g, std::ostringstream& o, bool& first) {
    if (g.getNumGeometries() > 1 || g.getGeometryTypeId() == geom::GEOS_GEOMETRYCOLLECTION ||
        g.getGeometryTypeId() == geom::GEOS_MULTILINESTRING || g.getGeometryTypeId() == geom::GEOS_MULTICURVE) {
        for (std::size_t i = 0; i < g.getNumGeometries(); i++) dumpGeom(*g.getGeometryN(i), o, first);
        return;
    }
    if (!first) o << ",";
    first = false;
    o << "{\"t\":\"" << g.getGeometryType() << "\"";
    if (g.getGeometryTypeId() == geom::GEOS_COMPOUNDCURVE) {
        const auto& cc = static_cast<const geom::Curve&>(g);
        o << ",\"parts\":[";
        for (std::size_t i = 0; i < cc.getNumCurves(); i++) {
            const auto* sc = cc.getCurveN(i);
            if (i) o << ",";
            o << "{\"t\":\"" << sc->getGeometryType() << "\",\"c\":[";
            const auto* cs = sc->getCoordinatesRO();
            for (std::size_t j = 0; j < cs->size(); j++) { if (j) o << ","; o << pt(cs->getAt<CoordinateXY>(j)); }
            o << "]}";
        }
        o << "]";
    } else {
        auto cs = g.getCoordinates();
        o << ",\"c\":[";
        for (std::size_t j = 0; j < cs->size(); j++) { if (j) o << ","; o << pt(cs->getAt<CoordinateXY>(j)); }
        o << "]";
    }
    o << "}";
}

static std::string handle(const std::string& line) {
    std::istringstream in(line);
    std::string op; in >> op;
    std::ostringstream o;
    if (op == "A") {
        auto a = rp(in), b = rp(in), c = rp(in); double dev = rd(in);
        auto s = seq3(a, b, c);
        CircularArc arc(*s, 0);
        auto r = seq3(c, b, a);
        CircularArc rev(*r, 0);
        o << "{\"op\":\"A\"," << arcProps(arc, dev) << ",\"rev\":{" << arcProps(rev, dev) << "}}";
    } else if (op == "AA") {
        auto a0 = rp(in), a1 = rp(in), a2 = rp(in), b0 = rp(in), b1 = rp(in), b2 = rp(in);
        auto sa = seq3(a0, a1, a2), sb = seq3(b0, b1, b2);
        CircularArc A(*sa, 0), B(*sb, 0);
        CircularArcIntersector ci1, ci2;
        ci1.intersects(A, B);
        ci2.intersects(B, A);
        auto cp = A.closestPoints(B);
        o << "{\"op\":\"AA\",\"ab\":" << ciResult(ci1) << ",\"ba\":" << ciResult(ci2)
          << ",\"dist\":" << num(A.distance(B)) << ",\"cp\":[" << pt(cp[0]) << "," << pt(cp[1]) << "]}";
    } else if (op == "AS") {
        auto a0 = rp(in), a1 = rp(in), a2 = rp(in), s0 = rp(in), s1 = rp(in);
        auto sa = seq3(a0, a1, a2);
        CoordinateSequence seg(2u, false, false);
        seg.setAt(s0, 0); seg.setAt(s1, 1);
        CoordinateSequence segr(2u, false, false);
        segr.setAt(s1, 0); segr.setAt(s0, 1);
        CircularArc A(*sa, 0);
        CircularArcIntersector ci1, ci2;
        ci1.intersects(A, seg, 0, 1, false);
        ci2.intersects(A, segr, 0, 1, false);
        auto cp = A.closestPoints(s0, s1);
        o << "{\"op\":\"AS\",\"fwd\":" << ciResult(ci1) << ",\"rev\":" << ciResult(ci2)
          << ",\"dist\":" << num(A.distance(s0, s1)) << ",\"cp\":[" << pt(cp[0]) << "," << pt(cp[1]) << "]}";
    } else if (op == "AP") {
        auto a0 = rp(in), a1 = rp(in), a2 = rp(in), p = rp(in);
        auto sa = seq3(a0, a1, a2);
        CircularArc A(*sa, 0);
        o << "{\"op\":\"AP\",\"dist\":" << num(A.distance(p)) << ",\"cp\":" << pt(A.closestPoint(p))
          << ",\"contains\":" << (A.containsPoint(p) ? "true" : "false") << "}";
    } else if (op == "SS") {
        auto p0 = rp(in), p1 = rp(in), q0 = rp(in), q1 = rp(in);
        CoordinateSequence P(2u, false, false), Q(2u, false, false);
        P.setAt(p0, 0); P.setAt(p1, 1); Q.setAt(q0, 0); Q.setAt(q1, 1);
        CircularArcIntersector ci;
        ci.intersects(P, 0, 1, Q, 0, 1);
        o << "{\"op\":\"SS\",\"r\":" << ciResult(ci) << "}";
    } else if (op == "G" || op == "N") {
        int preserve = 0;
        if (op == "N") in >> preserve;
        std::string wkt; std::getline(in, wkt);
        io::WKTReader reader;
        auto g = reader.read(wkt);
        o << "{\"op\":\"" << op << "\",\"type\":\"" << g->getGeometryType() << "\"";
        if (op == "G") {
            o << ",\"length\":" << num(g->getLength()) << ",\"area\":" << num(g->getArea());
            auto lin = g->getLinearized(algorithm::CurveToLineParams::maxDeviation(1e-3 * std::max(1e-300, g->getEnvelopeInternal()->getDiameter())));
            o << ",\"linLength\":" << num(lin->getLength()) << ",\"linArea\":" << num(lin->getArea())
              << ",\"linType\":\"" << lin->getGeometryType() << "\"";
        } else {
            noding::GeometryNoder noder(*g);
            // setPreserveCompoundCurves(true) asserts onlyFirstGeomEdges (GeometryNoder.cpp:349);
            // GeometrySplitter, its only in-tree caller, always sets both.
            noder.setOnlyFirstGeomEdges(preserve != 0);
            noder.setPreserveCompoundCurves(preserve != 0);
            auto res = noder.getNoded();
            o << ",\"out\":[";
            bool first = true;
            dumpGeom(*res, o, first);
            o << "]";
        }
        o << "}";
    } else {
        o << "{\"op\":\"?\"}";
    }
    return o.str();
}

int main() {
    std::ios::sync_with_stdio(false);
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        try {
            std::cout << handle(line) << "\n";
        } catch (const std::exception& e) {
            std::string msg = e.what();
            for (auto& ch : msg) if (ch == '"' || ch == '\\' || ch == '\n') ch = ' ';
            std::cout << "{\"op\":\"EXC\",\"msg\":\"" << msg << "\"}\n";
        }
        std::cout.flush();
    }
    return 0;
}
