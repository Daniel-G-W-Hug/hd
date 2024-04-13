#pragma once

// author: Daniel Hug, 2024

#include <array>
#include <cmath>    // abs, sqrt, acos
#include <compare>  // <=>
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
// Vec2d<T> definition
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct Vec2d {

    // ctor
    Vec2d<T>() = default;
    Vec2d<T>(T x_in, T y_in) : x(x_in), y(y_in) {}

    T x{};
    T y{};

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(const Vec2d<U>& rhs) const
    {
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_x = std::abs(rhs.x - x);
        auto abs_delta_y = std::abs(rhs.y - y);
        auto constexpr delta_eps =
            3.0 * std::max<std::common_type_t<T, U>>(std::numeric_limits<T>::epsilon(),
                                                     std::numeric_limits<U>::epsilon());
        if (abs_delta_x < delta_eps && abs_delta_y < delta_eps) return true;
        return false;
    }

    // unary minus (must be declared a friend otherwise doesn't work)
    friend inline Vec2d<T> operator-(const Vec2d<T>& v) { return Vec2d<T>(-v.x, -v.y); }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Vec2d<U>& v);
};

// for printing via iostream
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, const Vec2d<T>& v)
{
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> operations
////////////////////////////////////////////////////////////////////////////////

// adding vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator+(const Vec2d<T>& v1,
                                                           const Vec2d<U>& v2)
{
    return Vec2d<std::common_type_t<T, U>>(v1.x + v2.x, v1.y + v2.y);
}

// substracting vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator-(const Vec2d<T>& v1,
                                                           const Vec2d<U>& v2)
{
    return Vec2d<std::common_type_t<T, U>>(v1.x - v2.x, v1.y - v2.y);
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(const Vec2d<T>& v, U s)
{
    return Vec2d<std::common_type_t<T, U>>(v.x * s, v.y * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(T s, const Vec2d<U>& v)
{
    return Vec2d<std::common_type_t<T, U>>(v.x * s, v.y * s);
}

// devide a vector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator/(const Vec2d<T>& v, U s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return Vec2d<std::common_type_t<T, U>>(v.x * inv, v.y * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of vector
template <typename T> inline T sq_nrm(const Vec2d<T>& v) { return v.x * v.x + v.y * v.y; }

// return magnitude of vector
template <typename T> inline T nrm(const Vec2d<T>& v)
{
    return std::sqrt(v.x * v.x + v.y * v.y);
}

// return a vector normalized to norm(v) == 1.0
template <typename T> inline Vec2d<T> normalized(const Vec2d<T>& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return Vec2d<T>(v.x * inv, v.y * inv);
}

// return the multiplicative inverse of the vector
template <typename T> inline Vec2d<T> inv(const Vec2d<T>& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = 1.0 / sq_n; // for multiplication with inverse of squared norm
    return Vec2d<T>(v.x * inv, v.y * inv);
}

// return dot-product of two vectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(const Vec2d<T>& v1, const Vec2d<U>& v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(const Vec2d<T>& v1, const Vec2d<U>& v2)
{
    std::common_type_t<T, U> norm_prod = nrm(v1) * nrm(v2);
    if (norm_prod < std::numeric_limits<std::common_type_t<T, U>>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(norm_prod) + "\n");
    }
    return std::acos(dot(v1, v2) / norm_prod);
}

// wedge product (returns a bivector, which is the pseudoscalar in 2d)
// wedge(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> wdg(const Vec2d<T>& v1, const Vec2d<U>& v2)
{
    return v1.x * v2.y - v1.y * v2.x;
}

// projection of v1 onto v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto(const Vec2d<T>& v1,
                                                              const Vec2d<U>& v2)
{
    return dot(v1, v2) * Vec2d<std::common_type_t<T, U>>(inv(v2));
}

// projection of v1 onto v2 (v2 must already be normalized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto_n(const Vec2d<T>& v1,
                                                                const Vec2d<U>& v2)
{
    return dot(v1, v2) * v2;
}

// rejection of v1 from v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from(const Vec2d<T>& v1,
                                                             const Vec2d<U>& v2)
{
    // version using vector substraction
    // return v1 - project_onto(v1, v2);

    using ctype = std::common_type_t<T, U>;
    // version using geometric algebra wedge product manually computed
    // from "wdg(v1,v2)*inv(v2)"
    ctype w = wdg(v1, v2); // bivector with component e12
    ctype sq_n = sq_nrm(v2);
    if (sq_n < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    ctype w_sq_n_inv = w / sq_n;
    return Vec2d<ctype>(v2.y * w_sq_n_inv, -v2.x * w_sq_n_inv);
}

// rejection of v1 from v2 (v2 must already be normalized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from_n(const Vec2d<T>& v1,
                                                               const Vec2d<U>& v2)
{
    return v1 - dot(v1, v2) * v2;
}

} // namespace hd::ga


////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T> struct fmt::formatter<hd::ga::Vec2d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec2d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({}, {})", nested(v.x), nested(v.y));
    }
};

// Usage:
//
// std::vector<Vec2d<double>> vp1{{1.0, 1.0}, {1.5, 2.0}};
// Vec2d p{1.0, 2.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));
