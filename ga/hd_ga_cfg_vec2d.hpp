#pragma once

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
    bool operator==(const Vec2d<T>& rhs) const
    {
        // comparison using componentwise comparison
        // equality implies same magnitude and direction
        auto adx = std::abs(rhs.x - x);
        auto ady = std::abs(rhs.y - y);
        if (adx < std::numeric_limits<T>::epsilon() &&
            ady < std::numeric_limits<T>::epsilon())
            return true;
        return false;
    }

    // unary minus
    inline Vec2d<T> operator-(const Vec2d<T>& v) { return Vec2d<T>(-v.x, -v.y); }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Vec2d<U>& v);
};

// for printing via iostream
template <typename T> std::ostream& operator<<(std::ostream& os, const Vec2d<T>& v)
{
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> operations
////////////////////////////////////////////////////////////////////////////////

// adding vectors
template <typename T> inline Vec2d<T> operator+(const Vec2d<T>& v1, const Vec2d<T>& v2)
{
    return Vec2d<T>(v1.x + v2.x, v1.y + v2.y);
}

// substracting vectors
template <typename T> inline Vec2d<T> operator-(const Vec2d<T>& v1, const Vec2d<T>& v2)
{
    return Vec2d<T>(v1.x - v2.x, v1.y - v2.y);
}


// multiply a vector with a scalar
template <typename T> inline Vec2d<T> operator*(const Vec2d<T>& v, T s)
{
    return Vec2d<T>(v.x * s, v.y * s);
}

template <typename T> inline Vec2d<T> operator*(T s, const Vec2d<T>& v)
{
    return Vec2d<T>(v.x * s, v.y * s);
}

// devide a vector by a scalar
template <typename T> inline Vec2d<T> operator/(const Vec2d<T>& v, T s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    T inv = 1.0 / s; // for multiplicaton with inverse value
    return Vec2d<T>(v.x * inv, v.y * inv);
}

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of vector
template <typename T> inline T sq_norm(const Vec2d<T>& v)
{
    return v.x * v.x + v.y * v.y;
}

// return magnitude of vector
template <typename T> inline T norm(const Vec2d<T>& v)
{
    return std::sqrt(v.x * v.x + v.y * v.y);
}

// return a vector normalized to norm(v) == 1.0
template <typename T> inline Vec2d<T> normalized(const Vec2d<T>& v)
{
    T n = norm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return Vec2d<T>(v.x * inv, v.y * inv);
}

// return the multiplicative inverse of the vector
template <typename T> inline Vec2d<T> inverse(const Vec2d<T>& v)
{
    T sq_n = sq_norm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = 1.0 / sq_n; // for multiplication with inverse of squared norm
    return Vec2d<T>(v.x * inv, v.y * inv);
}

// return dot-product of two vectors
template <typename T> inline T dot(const Vec2d<T>& v1, const Vec2d<T>& v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T> inline T angle(const Vec2d<T>& v1, const Vec2d<T>& v2)
{
    T norm_prod = norm(v1) * norm(v2);
    if (norm_prod < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(norm_prod) + "\n");
    }
    return std::acos(dot(v1, v2) / norm_prod);
}

// wedge product (returns a bivector, which is the pseudoscalar in 2d)
// wedge(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T> inline T wedge(const Vec2d<T>& v1, const Vec2d<T>& v2)
{
    return v1.x * v2.y - v1.y * v2.x;
}

// projection of v1 onto v2
template <typename T> inline Vec2d<T> project_onto(const Vec2d<T>& v1, const Vec2d<T>& v2)
{
    return dot(v1, v2) * Vec2d<T>(inverse(v2));
}

// projection of v1 onto v2 (v2 must already be normalized to norm(v2) == 1)
template <typename T>
inline Vec2d<T> project_onto_n(const Vec2d<T>& v1, const Vec2d<T>& v2)
{
    return dot(v1, v2) * v2;
}

// rejection of v1 from v2
template <typename T> inline Vec2d<T> reject_from(const Vec2d<T>& v1, const Vec2d<T>& v2)
{
    // version using vector substraction
    // return v1 - project_onto(v1, v2);

    // version using geometric algebra wedge product manually computed
    // from "wedge(v1,v2)*inverse(v2)"
    T w = wedge(v1, v2); // bivector with component e12
    T sq_n = sq_norm(v2);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T w_sq_n_inv = w / sq_n;
    return Vec2d<T>(v2.y * w_sq_n_inv, -v2.x * w_sq_n_inv);
}

// rejection of v1 from v2 (v2 must already be normalized to norm(v2) == 1)
template <typename T>
inline Vec2d<T> reject_from_n(const Vec2d<T>& v1, const Vec2d<T>& v2)
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
