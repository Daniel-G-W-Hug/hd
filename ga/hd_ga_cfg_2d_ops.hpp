#pragma once

// author: Daniel Hug, 2024

#include <cmath>    // abs, sqrt, acos
#include <concepts> // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "hd_ga_cfg_value_t.hpp"

#include "hd_ga_cfg_mvec2d.hpp"
#include "hd_ga_cfg_vec2d.hpp"

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return dot-product of two vectors
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(const Vec2d<T>& v1, const Vec2d<U>& v2)
{
    // this implementation is only valid in an orthonormal basis
    return v1.x * v2.x + v1.y * v2.y;
}

// return squared magnitude of vector
template <typename T> inline T sq_nrm(const Vec2d<T>& v) { return dot(v, v); }

// return magnitude of vector
template <typename T> inline T nrm(const Vec2d<T>& v) { return std::sqrt(dot(v, v)); }

// return a vector unitized to nrm(v) == 1.0
template <typename T> inline Vec2d<T> unitized(const Vec2d<T>& v)
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
    T inv = 1.0 / sq_n; // inverse of squared norm for a vector
    return Vec2d<T>(v.x * inv, v.y * inv);
}

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(const Vec2d<T>& v1, const Vec2d<U>& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    return std::acos(dot(v1, v2) / nrm_prod);
}

// wedge product (returns a bivector, which is the pseudoscalar in 2d)
// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar2d<std::common_type_t<T, U>> wdg(const Vec2d<T>& v1, const Vec2d<U>& v2)
{
    return PScalar2d<std::common_type_t<T, U>>(v1.x * v2.y - v1.y * v2.x);
}

// projection of v1 onto v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto(const Vec2d<T>& v1,
                                                              const Vec2d<U>& v2)
{
    return dot(v1, v2) * inv(v2);
}

// projection of v1 onto v2 (v2 must already be unitized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto_unitized(const Vec2d<T>& v1,
                                                                       const Vec2d<U>& v2)
{
    // requires v2 to be unitized
    return dot(v1, v2) * v2;
}

// rejection of v1 from v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from(const Vec2d<T>& v1,
                                                             const Vec2d<U>& v2)
{
    using ctype = std::common_type_t<T, U>;
    // version using geometric algebra wedge product manually computed
    // from "wdg(v1,v2)*inv(v2)"
    PScalar2d<ctype> w = wdg(v1, v2); // bivector with component e12
    ctype sq_n = sq_nrm(v2);          //
    if (sq_n < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    ctype w_sq_n_inv = w / sq_n;
    return Vec2d<ctype>(v2.y * w_sq_n_inv, -v2.x * w_sq_n_inv);
}

// rejection of v1 from v2 (v2 must already be unitized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from_unitized(const Vec2d<T>& v1,
                                                                      const Vec2d<U>& v2)
{
    // requires v2 to be unitized

    using ctype = std::common_type_t<T, U>;
    // version using geometric algebra wedge product manually computed
    // from "wdg(v1,v2)*inv(v2)" + v2 being already it's own inverse
    PScalar2d<ctype> w = wdg(v1, v2); // bivector with component e12
    return Vec2d<ctype>(v2.y * w, -v2.x * w);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// geometric product ab for fully populated 2d multivector
// gpr() ... geometric product
// Expensive! - Don't use if you don't have to! (16x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab) (vector vector = scalar + bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MVec2d<T> const& v1,
                                                      MVec2d<U> const& v2)
{
    // geometric product with a fully populated 2d multivector
    T c0 = v1.c0 * v2.c0 + v1.c1 * v2.c1 + v1.c2 * v2.c2 - v1.c3 * v2.c3;
    T c1 = v1.c0 * v2.c1 + v1.c1 * v2.c0 - v1.c2 * v2.c3 + v1.c3 * v2.c2;
    T c2 = v1.c0 * v2.c2 + v1.c1 * v2.c3 + v1.c2 * v2.c0 - v1.c3 * v2.c1;
    T c3 = v1.c0 * v2.c3 + v1.c1 * v2.c2 - v1.c2 * v2.c1 + v1.c3 * v2.c0;
    return MVec2d<std::common_type_t<T, U>>(c0, c1, c2, c3);
}

// geometric product ab for two vectors (returns a multivector)
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab) (vector vector = scalar + bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(Vec2d<T> const& a,
                                                      Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(Scalar<ctype>(dot(a, b)), wdg(a, b));
}

// geometric product AB of a trivector A multiplied from the left
// to the multivector B
// gpr(trivector, multivector) => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(PScalar2d<T> A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MVec2d<ctype>(-B.c3, B.c2, -B.c1, B.c0);
}

// geometric product Ab of a trivector A multiplied from the left
// to the vector b
// gpr(trivector, vector) => returns a vector
// this multiplication turns the vector by -90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> gpr(PScalar2d<T> A, Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * Vec2d<ctype>(b.y, -b.x);
}

// geometric product AB of a trivector B multiplied from the right
// to the multivector A
// gpr(multivector, trivector) => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MVec2d<T> const& A, PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-A.c3, -A.c2, A.c1, A.c0);
}

// geometric product aB of a trivector B multiplied from the right
// to the vector a
// gpr(vector, trivector) => returns a vector
// this multiplication turns the vector by +90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> gpr(Vec2d<T> const& a, PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-a.y, a.x) * ctype(B);
}

// geometric product AB of two bivectors
// gpr(bivector, bivector) => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> gpr(PScalar2d<T> A, PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * ctype(B); // bivectors square to -1
}

} // namespace hd::ga
