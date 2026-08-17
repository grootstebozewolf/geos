//
// Test Suite for geos::algorithm::RocqNative (Phase 5 libntsrocq FFI).
// Skips kernel checks when the shared library is not loadable.

#include <tut/tut.hpp>
#include <geos/algorithm/RocqNative.h>

namespace tut {
struct test_rocqnative_data {};
typedef test_group<test_rocqnative_data> group;
typedef group::object object;

group test_rocqnative_group("geos::algorithm::RocqNative");

// testUnavailableIsSafe
template<>
template<>
void object::test<1>()
{
    // Must not throw just because the native library is absent.
    (void) ntsrocq::RocqNative::isAvailable();
}

// testOrientFilteredCCW
template<>
template<>
void object::test<2>()
{
    if (!ntsrocq::RocqNative::isAvailable()) {
        return;
    }
    auto s = ntsrocq::RocqNative::orientSignFiltered(0, 0, 1, 0, 0, 1);
    ensure_equals(static_cast<int>(s), static_cast<int>(ntsrocq::OrientSign::Pos));
}

// testInCircleOutside
template<>
template<>
void object::test<3>()
{
    if (!ntsrocq::RocqNative::isAvailable()) {
        return;
    }
    double det = ntsrocq::RocqNative::inCircle(0, 0, 2, 0, 1, 1, 1, -0.5);
    ensure(det < 0);
}

} // namespace tut
