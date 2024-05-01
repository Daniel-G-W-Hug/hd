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


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> definition (used for implementation of algebra<3,0,0>)
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct Vec3d {

    // assumes a right-handed orthonormal vector basis {e1, e2, e3}
    // using components {x, y, z}, such that for each vector v:
    // v = x * e1 + y * e2 + z * e3

    // ctors
    Vec3d() = default;

    Vec3d(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {}

    template <typename U>
        requires(std::floating_point<U>)
    Vec3d(Vec3d<U> const& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    T x{};
    T y{};
    T z{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(Vec3d<U> const& rhs) const
    {
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_x = std::abs(rhs.x - x);
        auto abs_delta_y = std::abs(rhs.y - y);
        auto abs_delta_z = std::abs(rhs.z - z);
        auto constexpr delta_eps =
            std::common_type_t<T, U>(5.0) *
            std::max<std::common_type_t<T, U>>(std::numeric_limits<T>::epsilon(),
                                               std::numeric_limits<U>::epsilon());
        if (abs_delta_x < delta_eps && abs_delta_y < delta_eps && abs_delta_z < delta_eps)
            return true;
        return false;
    }

    // unary minus (must be declared a friend otherwise doesn't work)
    friend inline constexpr Vec3d<T> operator-(Vec3d<T> const& v)
    {
        return Vec3d<T>(-v.x, -v.y, -v.z);
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, Vec3d<U> const& v);
};

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// adding vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator+(Vec3d<T> const& v1,
                                                           Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

// substracting vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator-(Vec3d<T> const& v1,
                                                           Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v, U s)
{
    return Vec3d<std::common_type_t<T, U>>(v.x * s, v.y * s, v.z * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(T s, Vec3d<U> const& v)
{
    return Vec3d<std::common_type_t<T, U>>(v.x * s, v.y * s, v.z * s);
}

// devide a vector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator/(Vec3d<T> const& v, U s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return Vec3d<std::common_type_t<T, U>>(v.x * inv, v.y * inv, v.z * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return the dot product of two vectors (= a scalar)
// coordinate free definition: dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    // this implementation is only valid in an orthonormal basis
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// return squared magnitude of vector
template <typename T> inline T sq_nrm(Vec3d<T> const& v) { return dot(v, v); }

// return magnitude of vector
template <typename T> inline T nrm(Vec3d<T> const& v) { return std::sqrt(dot(v, v)); }

// return a vector unitized to nrm(v) == 1.0
template <typename T> inline Vec3d<T> unitized(Vec3d<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return Vec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the multiplicative inverse of the vector
template <typename T> inline Vec3d<T> inv(Vec3d<T> const& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = 1.0 / sq_n; // inverse of squared norm for a vector
    return Vec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    std::common_type_t<T, U> nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<std::common_type_t<T, U>>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    return std::acos(dot(v1, v2) / nrm_prod);
}

// cross-product between two vectors (returns a vector in 3d)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> cross(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Vec3d<T> const& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
}

} // namespace hd::ga


////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T> struct fmt::formatter<hd::ga::Vec3d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec3d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{})", nested(v.x), nested(v.y),
                              nested(v.z));
    }
};

// Usage:
//
// std::vector<Vec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
// Vec3d p{1.0, 2.0, 3.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));
