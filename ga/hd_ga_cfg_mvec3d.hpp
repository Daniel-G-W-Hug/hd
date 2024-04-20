#pragma once

// author: Daniel Hug, 2024

#include <algorithm>
#include <cmath>    // abs, sqrt, acos
#include <concepts> // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>


#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "hd_ga_cfg_value_t.hpp"
#include "hd_ga_cfg_vec2d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> definition of a multivector
////////////////////////////////////////////////////////////////////////////////

template <typename T = value_t>
    requires(std::floating_point<T>)
struct MVec3d {

    // ctors

    // (all grades = 0)
    MVec3d() = default;

    // assign all components directly
    MVec3d(T s, T x, T y, T z, T yz, T zx, T xy, T ps) :
        c0(s), c1(x), c2(y), c3(z), c4(yz), c5(zx), c6(xy), c7(ps)
    {
    }

    // assign a scalar part exclusively (other grades = 0)
    MVec3d(Scalar_t<T> s) : c0(s) {}

    // assign a vector part exclusively (other grades = 0)
    MVec3d(Vec3d<T> const& v) : c1(v.x), c2(v.y), c3(v.z) {}

    // assign a bivector part exclusively (other grades = 0)
    MVec3d(BiVec3d<T> const& v) : c4(v.x), c5(v.y), c6(v.z) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to computing the full geometric product)
    MVec3d(Scalar_t<T> s, BiVec3d<T> const& v) : c0(s), c4(v.x), c5(v.y), c6(v.z) {}

    // assign a geometric product resulting from a product of a vector and a bivector
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to computing the full geometric product)
    MVec3d(Vec3d<T> const& v, PScalar3d_t<T> ps) : c1(v.x), c2(v.y), c3(v.z), c7(ps) {}

    // this constructor must be explicitly deleted,
    // otherwise a BiVec3d<T> could implicitly convert to Vec3d<T>
    // reason: Vec3d<T> is a base class of BiVec3d<T>
    MVec3d(BiVec3d<T> const& v, PScalar3d_t<T> ps) = delete;

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec3d(PScalar3d_t<T> ps) : c7(ps) {}

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MVec3d<T>(MVec3d<U> const& v) :
        c0(v.c0), c1(v.c1), c2(v.c2), c3(v.c3), c4(v.c4), c5(v.c5), c6(v.c6), c7(v.c7)
    {
    }


    T c0{}; // scalar
    T c1{}; // vector 3d, 1st component   (x)  - maps to basis bivector  e1
    T c2{}; // vector 3d, 2nd component   (y)  - maps to basis bivector  e2
    T c3{}; // vector 3d, 3rd component   (z)  - maps to basis bivector  e3
    T c4{}; // bivector 3d, 1st component (yz) - maps to basis bivector  e2^e3
    T c5{}; // bivector 3d, 2nd component (zx) - maps to basis bivector  e3^e1
    T c6{}; // bivector 3d, 3rd component (xy) - maps to basis bivector  e1^e2
    T c7{}; // trivector 3d = 3d pseudoscalar  - maps to basis trivector e1^e2^e3

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(const MVec3d<U>& rhs) const
    {
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto abs_delta_c2 = std::abs(rhs.c2 - c2);
        auto abs_delta_c3 = std::abs(rhs.c3 - c3);
        auto abs_delta_c4 = std::abs(rhs.c4 - c4);
        auto abs_delta_c5 = std::abs(rhs.c5 - c5);
        auto abs_delta_c6 = std::abs(rhs.c6 - c6);
        auto abs_delta_c7 = std::abs(rhs.c7 - c7);
        auto constexpr delta_eps =
            std::common_type_t<T, U>(5.0) *
            std::max<std::common_type_t<T, U>>(std::numeric_limits<T>::epsilon(),
                                               std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps &&
            abs_delta_c2 < delta_eps && abs_delta_c3 < delta_eps &&
            abs_delta_c4 < delta_eps && abs_delta_c5 < delta_eps &&
            abs_delta_c6 < delta_eps && abs_delta_c7 < delta_eps)
            return true;
        return false;
    }

    // unary minus (must be declared a friend otherwise doesn't work)
    friend inline MVec3d<T> operator-(const MVec3d<T>& v)
    {
        return MVec3d<T>(-v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, -v.c7);
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const MVec3d<U>& v);
};

// for printing via iostream
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, const MVec3d<T>& v)
{
    os << "(" << v.c0 << ", " << v.c1 << ", " << v.c2 << ", " << v.c3 << ", " << v.c4
       << ", " << v.c5 << ", " << v.c6 << ", " << v.c7 << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// add multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator+(const MVec3d<T>& v1,
                                                            const MVec3d<U>& v2)
{
    return MVec3d<std::common_type_t<T, U>>(v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2,
                                            v1.c3 + v2.c3, v1.c4 + v2.c4, v1.c5 + v2.c5,
                                            v1.c6 + v2.c6, v1.c7 + v2.c7);
}

// substract multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator-(const MVec3d<T>& v1,
                                                            const MVec3d<U>& v2)
{
    return MVec3d<std::common_type_t<T, U>>(v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2,
                                            v1.c3 - v2.c3, v1.c4 - v2.c4, v1.c5 - v2.c5,
                                            v1.c6 - v2.c6, v1.c7 - v2.c7);
}


// multiply a multivector with a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(const MVec3d<T>& v, U s)
{
    return MVec3d<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s,
                                            v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(T s, const MVec3d<U>& v)
{
    return MVec3d<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s,
                                            v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s);
}

// devide a multivector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator/(const MVec3d<T>& v, U s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return MVec3d<std::common_type_t<T, U>>(v.c0 * inv, v.c1 * inv, v.c2 * inv,
                                            v.c3 * inv, v.c4 * inv, v.c5 * inv,
                                            v.c6 * inv, v.c7 * inv);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector (= pseudoscalar in 3d)

template <typename T> inline constexpr Scalar_t<T> gr0(const MVec3d<T>& v)
{
    return Scalar_t<T>(v.c0);
}

template <typename T> inline constexpr Vec3d<T> gr1(const MVec3d<T>& v)
{
    return Vec3d<T>(v.c1, v.c2, v.c3);
}
template <typename T> inline constexpr BiVec3d<T> gr2(const MVec3d<T>& v)
{
    return BiVec3d<T>(v.c4, v.c5, v.c6);
}

template <typename T> inline constexpr PScalar3d_t<T> gr3(const MVec3d<T>& v)
{
    return PScalar3d_t<T>(v.c7);
}


// geometric product AB for fully populated 3d multivector
// gpr() ... geometric product
// Expensive! - Don't use if you don't have to! (64x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector vector = scalar + bivector)
// Ab = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector vector = vector + trivector)
// aB = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)  (bivector vector = vector + trivector)
// => see overloaded versions of gpr
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> gpr(const MVec3d<T>& v1,
                                                      const MVec3d<U>& v2)
{
    // geometric product with a fully populated 3d multivector
    T c0 = v1.c0 * v2.c0 + v1.c1 * v2.c1 + v1.c2 * v2.c2 + v1.c3 * v2.c3 - v1.c4 * v2.c4 -
           v1.c5 * v2.c5 - v1.c6 * v2.c6 - v1.c7 * v2.c7;
    T c1 = v1.c0 * v2.c1 + v1.c1 * v2.c0 - v1.c2 * v2.c6 + v1.c3 * v2.c5 - v1.c4 * v2.c7 -
           v1.c5 * v2.c3 + v1.c6 * v2.c2 - v1.c7 * v2.c4;
    T c2 = v1.c0 * v2.c2 + v1.c1 * v2.c6 + v1.c2 * v2.c0 - v1.c3 * v2.c4 + v1.c4 * v2.c3 -
           v1.c5 * v2.c7 - v1.c6 * v2.c1 - v1.c7 * v2.c5;
    T c3 = v1.c0 * v2.c3 - v1.c1 * v2.c5 + v1.c2 * v2.c4 + v1.c3 * v2.c0 - v1.c4 * v2.c2 +
           v1.c5 * v2.c1 - v1.c6 * v2.c7 - v1.c7 * v2.c6;
    T c4 = v1.c0 * v2.c4 + v1.c1 * v2.c7 + v1.c2 * v2.c3 - v1.c3 * v2.c2 + v1.c4 * v2.c0 -
           v1.c5 * v2.c6 + v1.c6 * v2.c5 + v1.c7 * v2.c1;
    T c5 = v1.c0 * v2.c5 - v1.c1 * v2.c3 + v1.c2 * v2.c7 + v1.c3 * v2.c1 + v1.c4 * v2.c6 +
           v1.c5 * v2.c0 - v1.c6 * v2.c4 + v1.c7 * v2.c2;
    T c6 = v1.c0 * v2.c6 + v1.c1 * v2.c2 - v1.c2 * v2.c1 + v1.c3 * v2.c7 - v1.c4 * v2.c5 +
           v1.c5 * v2.c4 + v1.c6 * v2.c0 + v1.c7 * v2.c3;
    T c7 = v1.c0 * v2.c7 + v1.c1 * v2.c4 + v1.c2 * v2.c5 + v1.c3 * v2.c6 + v1.c4 * v2.c1 +
           v1.c5 * v2.c2 + v1.c6 * v2.c3 + v1.c7 * v2.c0;
    return MVec3d<std::common_type_t<T, U>>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product ab between two vectors (returns a multivector)
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector vector = scalar + bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> gpr(const Vec3d<T>& a,
                                                      const Vec3d<U>& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(Scalar_t<ctype>(dot(a, b)), wdg(a, b));
}

// geometric product Ab for a bivector and a vector (returns a multivector)
// Ab = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector vector = vector + trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> gpr(const BiVec3d<T>& A,
                                                      const Vec3d<U>& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(dot(A, b), wdg(A, b));
}

// geometric product Ab for a bivector and a vector (returns a multivector)
// aB = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)  (bivector vector = vector + trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> gpr(const Vec3d<T>& a,
                                                      const BiVec3d<U>& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(dot(a, B), wdg(a, B));
}

} // namespace hd::ga


// ////////////////////////////////////////////////////////////////////////////////
// // printing support via fmt library
// ////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MVec3d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec3d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({}, {}, {}, {}, {}, {}, {}, {})", nested(v.c0),
                              nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
                              nested(v.c5), nested(v.c6), nested(v.c7));
    }
};
