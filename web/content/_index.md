---
title: "GEOS"
date: 2021-10-04T14:13:20-07:00
draft: false
geekdocNav: true
---

GEOS is a C/C++ library for [computational geometry](https://en.wikipedia.org/wiki/Computational_geometry) with a focus on algorithms used in
[geographic information systems](https://en.wikipedia.org/wiki/Geographic_information_system) (GIS) software.
It implements the [OGC](https://www.ogc.org/standards/sfa)
[*Simple Features*](https://en.wikipedia.org/wiki/Simple_Features) geometry model and provides all the spatial functions in that standard as well as many others.
GEOS is a core dependency of [PostGIS](https://postgis.net), [QGIS](https://qgis.org), [GDAL](https://gdal.org), [Shapely](https://shapely.readthedocs.io/en/stable/project.html) and [many others]({{< ref "usage/bindings" >}}).

## Capabilities

Spatial Model and Functions

* **Geometry Model**: Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection
* **Predicates**: Intersects, Touches, Disjoint, Crosses, Within, Contains, Overlaps, Equals, Covers
* **Operations**: Union, Distance, Intersection, Symmetric Difference, Convex Hull, Envelope, Buffer, Simplify, Polygon Assembly, Valid, Area, Length,
* **Prepared geometry** (using internal spatial indexes)
* **Spatial Indexes**: STR (Sort-Tile-Recursive) packed R-tree spatial index
* **Input/Output**: OGC Well Known Text (WKT) and Well Known Binary (WKB) readers and writers.

## API Features

* [C API]({{< ref "usage/c_api" >}}) (provides long-term API and ABI stability)
* [C++ API]({{< ref "usage/cpp_api" >}}) (will likely change across versions)
* Thread safety (using the [reentrant C API]({{<ref "usage/c_api/#reentrantthreadsafe-api" >}}))

## License

GEOS is [open source software](https://opensource.com/resources/what-open-source) available under the terms of [GNU Lesser General Public License](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html) (LGPL).


## Governance

The GEOS project is run by a
[Project Steering Committee]({{< ref "psc" >}}) made up of developers
and contributors to the project and is a project of
[OSGeo](https://www.osgeo.org/projects/geos/). Major decisions are made through the
[RFC process]({{< ref "rfcs" >}}), where proposals are documented and then voted on by the committee.
We strive to keep our day-to-day interactions collegial and welcoming within a basic [Code of Conduct]({{< ref "coc" >}}).

GEOS is developed as part of a family of closely related, direct ports of the
same computational geometry code base, and coordinates with its sibling projects
([JTS](https://github.com/locationtech/jts/) and [NTS](https://github.com/NetTopologySuite/NetTopologySuite))
on algorithms, behaviour and test data.

## The JTS / GEOS / NTS Ecosystem

GEOS started as a direct port to C++ of the
[JTS Topology Suite](https://github.com/locationtech/jts/) (JTS),
and remains tightly bound to that project. Most core algorithms have
been developed in JTS and ported to GEOS when complete.

[NetTopologySuite](https://github.com/NetTopologySuite/NetTopologySuite) (NTS)
is a direct port of JTS to .NET and is part of the same ecosystem:

| Project | Language | Role |
|---------|----------|------|
| [JTS](https://github.com/locationtech/jts/) | Java | Reference implementation, where most new algorithms are developed |
| GEOS | C/C++ | Direct port of JTS to C++, with a stable C API used by PostGIS, QGIS, GDAL, Shapely and others |
| [NTS](https://github.com/NetTopologySuite/NetTopologySuite) | C# / .NET | Direct port of JTS to .NET, used by the .NET geospatial stack |

The three projects follow the same design and share algorithms, so that fixes and
improvements made in one can be carried over to the others.

The projects share testing data.  This helps to ascertain when failures
are caused by differences in implementation (one port fails and the others do not)
and when they are caused by algorithm (all libraries fail).

## Citation

To cite GEOS in publications use:

> GEOS contributors (2025). GEOS computational geometry library. Open Source Geospatial Foundation. URL https://libgeos.org/.

A BibTeX entry for LaTeX users is
```
@Manual{,
  title = {{GEOS} computational geometry library},
  author = {{GEOS contributors}},
  organization = {Open Source Geospatial Foundation},
  year = {2025},
  url = {https://libgeos.org/},
  doi = {10.5281/zenodo.11396894}
}
```
