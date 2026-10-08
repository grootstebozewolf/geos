#!/usr/bin/env bash
# Build the GEOS probe against an existing GEOS build tree (default: ../../build).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BUILD="${GEOS_BUILD:-$ROOT/build}"
mkdir -p "$HERE/bin"
g++ -std=c++17 -O2 -g -I"$ROOT/include" -I"$BUILD/include" "$HERE/geos_probe.cpp" \
    -L"$BUILD/lib" -Wl,-rpath,"$BUILD/lib" -lgeos -o "$HERE/bin/geos_probe"
echo "built $HERE/bin/geos_probe"
