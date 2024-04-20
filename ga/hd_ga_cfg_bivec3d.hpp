#pragma once

// author: Daniel Hug, 2024

#include <cmath>    // abs
#include <concepts> // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "hd_ga_cfg_value_t.hpp"
#include "hd_ga_cfg_vec3d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// BiVec3d<T> definition (used for implementation of algebra<3,0,0>)
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct BiVec3d : public Vec3d<T> {

    // assumes a right-handed orthonormal vector basis {e1, e2, e3}
    // using components {x, y, z}, such that for each vector v:
    // v = x * e1 + y * e2 + z * e3
    //
    // and for each bivector bv:
    // bv = yz * e2^e3 + zx * e3^e1 + xy * e1^e2
    //    =  x * e2^e3 +  y * e3^e1 +  z * e1^e2
    // (same names (x, y, z), but same semantic as one live above (yz, zx, xy)

    // this is a mapping of the components as aliases
    // such that vector components x, y and z correspond to the
    // normals of the corresponding plane elements represented by
    // bivector components yz, zx and xy
    // i.e. they can by converted to each other by a duality transformation

    // if we would introduce aliases via references,
    // we would double the memory needed for each bivector
    // => So we don't! And just go on with using x, y, z
    //
    // T& yz = Vec3d<T>::x; // maps to basis bivector e2^e3
    // T& zx = Vec3d<T>::y; // maps to basis bivector e3^e1
    // T& xy = Vec3d<T>::z; // maps to basis bivector e1^e2

    // => everything is directly re-used from Vec3d<T> w/o modification.

    BiVec3d() = default;

    BiVec3d(T x_in, T y_in, T z_in) : Vec3d<T>(x_in, y_in, z_in) {}

    template <typename U>
        requires(std::floating_point<U>)
    BiVec3d(BiVec3d<U> const& v) : Vec3d<U>(reinterpret_cast<Vec3d<U> const&>(v))
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// BiVec3d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// adding vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator+(const BiVec3d<T>& v1,
                                                             const BiVec3d<U>& v2)
{
    return BiVec3d<std::common_type_t<T, U>>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

// substracting vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator-(const BiVec3d<T>& v1,
                                                             const BiVec3d<U>& v2)
{
    return BiVec3d<std::common_type_t<T, U>>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(const BiVec3d<T>& v, U s)
{
    return BiVec3d<std::common_type_t<T, U>>(v.x * s, v.y * s, v.z * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(T s, const BiVec3d<U>& v)
{
    return BiVec3d<std::common_type_t<T, U>>(v.x * s, v.y * s, v.z * s);
}

// devide a vector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator/(const BiVec3d<T>& v, U s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return BiVec3d<std::common_type_t<T, U>>(v.x * inv, v.y * inv, v.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// BiVec3d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, const BiVec3d<T>& v)
{
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

} // namespace hd::ga


////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::BiVec3d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::BiVec3d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({}, {}, {})", nested(v.x), nested(v.y),
                              nested(v.z));
    }
};

// Usage:
//
// std::vector<BiVec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
// BiVec3d p{1.0, 2.0, 3.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));
