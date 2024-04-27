#pragma once

// author: Daniel Hug, 2024

#include <cmath>    // abs, sqrt, acos
#include <concepts> // std::floating_point<T>
#include <iostream>
#include <limits>
#include <numbers> // math constants like pi
#include <stdexcept>
#include <string>

#include "hd_ga_cfg_value_t.hpp"

#include "hd_ga_cfg_vec2d.hpp"

#include "hd_ga_cfg_mcplx2d.hpp"
#include "hd_ga_cfg_mvec2d.hpp"


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
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MVec2d<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // geometric product with a fully populated 2d multivector
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c3;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product ab for two vectors (returns a multivector of the even subalgebra)
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
// => vector vector = scalar + bivector
//
// HINT: if a full multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> gpr(Vec2d<T> const& a,
                                                       Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MCplx2d<ctype>(Scalar<ctype>(dot(a, b)), wdg(a, b));
}

// define geometric multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                             Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product AB of a bivector A multiplied from the left
// to the multivector B
// gpr(bivector, multivector) => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(PScalar2d<T> A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MVec2d<ctype>(-B.c3, B.c2, -B.c1, B.c0);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> A,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product AB of a bivector A multiplied from the left
// to the multivector from the even subalgebra B (MCplx2d)
// gpr(bivector, even grade multivector) => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> gpr(PScalar2d<T> A,
                                                       MCplx2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MCplx2d<ctype>(-B.c1, B.c0);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator*(PScalar2d<T> A,
                                                             MCplx2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product Ab of a bivector A multiplied from the left
// to the vector b
// gpr(bivector, vector) => returns a vector
// this multiplication turns the vector by -90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> gpr(PScalar2d<T> A, Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * Vec2d<ctype>(b.y, -b.x);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> A,
                                                           Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, b);
}

// geometric product AB of a bivector B multiplied from the right
// to the multivector A
// gpr(multivector, bivector) => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MVec2d<T> const& A, PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-A.c3, -A.c2, A.c1, A.c0) * ctype(B);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product AB of a bivector B multiplied from the right
// to the even grade multivector A (MCplx2d)
// gpr(even grade multivector, bivector) => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> gpr(MCplx2d<T> const& A,
                                                       PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MCplx2d<ctype>(-A.c1, A.c0) * ctype(B);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator*(MCplx2d<T> const& A,
                                                             PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product aB of a bivector B multiplied from the right
// to the vector a
// gpr(vector, bivector) => returns a vector
// this multiplication turns the vector by +90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> gpr(Vec2d<T> const& a, PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-a.y, a.x) * ctype(B);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                           PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, B);
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

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> operator*(PScalar2d<T> A, PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

////////////////////////////////////////////////////////////////////////////////
// MCplx2d<T> operations
////////////////////////////////////////////////////////////////////////////////

// geometric product aB for a full 2d multivector B with a multivector from
// the even subalgebra a (only gr0 and gr2 components, i.e. a MCplx2d)
// => product with MCplx2d from the left
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MCplx2d<T> const& a,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(a.c0 * B.c0 - a.c1 * B.c3, a.c0 * B.c1 + a.c1 * B.c2,
                         a.c0 * B.c2 - a.c1 * B.c1, a.c0 * B.c3 + a.c1 * B.c0);
}

// define geometric multiplication with operator*(a,B) as an alias for gpr(a,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MCplx2d<T> const& a,
                                                            MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, B);
}

// geometric product ab for a 2d vector b with a multivector from
// the even subalgebra a (only gr0 and gr2 components, i.e. a MCplx2d)
// => product with MCplx2d from the left
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> gpr(MCplx2d<T> const& a,
                                                     Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(a.c0 * b.x + a.c1 * b.y, a.c0 * b.y - a.c1 * b.x);
}

// define geometric multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(MCplx2d<T> const& a,
                                                           Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product Ab for a full 2d multivector A with a multivector from
// the even subalgebra b (only gr0 and gr2 components, i.e. a MCplx2d)
// => product with MCplx2d from the right
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MVec2d<T> const& A,
                                                      MCplx2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c0 * b.c0 - A.c3 * b.c1, A.c1 * b.c0 - A.c2 * b.c1,
                         A.c1 * b.c1 + A.c2 * b.c0, A.c0 * b.c1 + A.c3 * b.c0);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            MCplx2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, b);
}

// geometric product aB for a full 2d vector a with a multivector from
// the even subalgebra b (only gr0 and gr2 components, i.e. a MCplx2d)
// => product with MCplx2d from the right
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> gpr(Vec2d<T> const& a,
                                                     MCplx2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(a.x * b.c0 - a.y * b.c1, a.x * b.c1 + a.y * b.c0);
}

// define geometric multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                           MCplx2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product ab for two multivectors from the even subalgebra (2d case)
// a  = gr0(a)  + gr2(a)  (scalar + bivector)
// b  = gr0(b)  + gr2(b)  (scalar + bivector)
// ab = gr0(ab) + gr2(ab) (scalar + bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> gpr(MCplx2d<T> const& a,
                                                       MCplx2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MCplx2d<ctype>(a.c0 * b.c0 - a.c1 * b.c1, a.c0 * b.c1 + a.c1 * b.c0);
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator*(MCplx2d<T> const& a,
                                                             MCplx2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// complex exponential function for setup of complex numbers
// as geometric multivector with a scalar and a bivector part
//
// r = 1 is the vector length of the complex number in polar form
// theta is the bivector angle (i.e. a multiple of the bivector I_2d)
// such that uv = r exp(I_2d theta) = a + I_2d b
// with r = |u| |v| = sqrt(a^2 + b^2) = 1
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MCplx2d<T> exp(PScalar2d<T> theta)
{
    return MCplx2d<T>(Scalar<T>(std::cos(theta)), PScalar2d<T>(std::sin(theta)));
}


// return squared magnitude of complex number
template <typename T> inline T sq_nrm(MCplx2d<T> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1;
}

// return magnitude of complex number
template <typename T> inline T nrm(MCplx2d<T> const& v) { return std::sqrt(sq_nrm(v)); }

// return conjugate complex of complex number,
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline MCplx2d<T> rev(MCplx2d<T> const& v)
{
    return MCplx2d<T>(v.c0, -v.c1);
}

// return a complex unitized to nrm(v) == 1.0
template <typename T> inline MCplx2d<T> unitized(MCplx2d<T> const& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return MCplx2d<T>(v.c0 * inv, v.c1 * inv);
}

// return the multiplicative inverse of the complex (inv(z) = 1/sq_nrm(z)*rev(z))
// with rev(z) being the complex conjugate
template <typename T> inline MCplx2d<T> inv(MCplx2d<T> const& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("complex norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = 1.0 / sq_n; // inverse of squared norm for a vector
    return inv * rev(v);
}

// return the angle of the complex number w.r.t. the real axis
// range of angle: -pi <= angle <= pi
template <typename T>
    requires(std::floating_point<T>)
inline T angle(MCplx2d<T> const& v)
{
    using std::numbers::pi;
    if (v.c0 > 0.0) {
        // quadrant I & IV
        return std::atan(v.c1 / v.c0);
    }
    if (v.c0 < 0.0 && v.c1 >= 0.0) {
        // quadrant II
        return std::atan(v.c1 / v.c0) + pi;
    }
    if (v.c0 < 0.0 && v.c1 < 0.0) {
        // quadrant III
        return std::atan(v.c1 / v.c0) - pi;
    }
    if (v.c0 == 0.0) {
        // on y-axis
        if (v.c1 > 0.0) return pi / 2.0;
        if (v.c1 < 0.0) return -pi / 2.0;
    }
    return 0.0; // zero as input => define 0 as corresponding angle
}

} // namespace hd::ga
