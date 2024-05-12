#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::clamp
#include <cmath>     // std::abs, std::sin, std::cos
#include <concepts>  // std::floating_point<T>
#include <iostream>
#include <limits>
#include <numbers> // math constants like pi
#include <stdexcept>
#include <string>

#include "hd_ga_cfg_value_t.hpp"

#include "hd_ga_cfg_vec3d.hpp"

#include "hd_ga_cfg_bivec3d.hpp"

#include "hd_ga_cfg_mvec3d.hpp"
#include "hd_ga_cfg_mvec3d_e.hpp"
#include "hd_ga_cfg_mvec3d_u.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// PScalar3d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of the pseudoscalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr T sq_nrm(PScalar3d<T> const& ps)
{
    return T(ps) * T(ps);
}

// return magnitude of the pseudoscalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr T nrm(PScalar3d<T> const& ps)
{
    return std::abs(T(ps));
}

// return inverse of the pseudoscalar (A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
// k is the dimension of the space of the pseudoscalar formed by k orthogonal vectors
template <typename T>
    requires(std::floating_point<T>)
inline constexpr T inv(PScalar3d<T> const& ps)
{
    return -T(ps) / sq_nrm(ps);
}

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> & BiVec3d<T> mixed geometric operations
////////////////////////////////////////////////////////////////////////////////

// return the dot product of a bivector and a vector (= a vector)
// dot(A,b) = gr0( gpr(A,b) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> dot(BiVec3d<T> const& v1, Vec3d<U> const& v2)
{
    // this implementation is only valid in an orthonormal basis
    return Vec3d<std::common_type_t<T, U>>(
        v1.z * v2.y - v1.y * v2.z, v1.x * v2.z - v1.z * v2.x, v1.y * v2.x - v1.x * v2.y);
}

// return the dot product of a vector and a bivector (= a vector)
// dot(a,B) = gr1( gpr(a,B) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> dot(Vec3d<T> const& v1, BiVec3d<U> const& v2)
{
    // this implementation is only valid in an orthonormal basis
    return Vec3d<std::common_type_t<T, U>>(
        v1.z * v2.y - v1.y * v2.z, v1.x * v2.z - v1.z * v2.x, v1.y * v2.x - v1.x * v2.y);
}


// return commutator product cmt(A,B) of two bivectors A and B (= a bivector)
// cmt(A,B) = 0.5*(AB-BA) = gr2( gpr(A,B) )
// the commutator product is antisymmetric, i.e. it is zero when a bivector is
// multiplied by itself, i.e. in that case only the dot product remains
// as the symmetric part
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec3d<std::common_type_t<T, U>> cmt(BiVec3d<T> const& A, BiVec3d<U> const& B)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(A.z * B.y - A.y * B.z, A.x * B.z - A.z * B.x,
                          A.y * B.x - A.x * B.y);
}

// return the angle between of a vector and a bivector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec3d<T> const& v1, BiVec3d<U> const& v2)
{
    std::common_type_t<T, U> nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<std::common_type_t<T, U>>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(dot(v1, v2) / nrm_prod, -1.0, 1.0));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(BiVec3d<T> const& v1, Vec3d<U> const& v2)
{
    std::common_type_t<T, U> nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<std::common_type_t<T, U>>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(dot(v1, v2) / nrm_prod, -1.0, 1.0));
}

// wedge product between two vectors (returns a bivector in 3d)
// coordinate-free definition: wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    return BiVec3d<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// wedge product between a vector and a bivector (returns a trivector in 3d)
// wdg(a,B) = gr3( gpr(a,B) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& a, BiVec3d<U> const& B)
{
    return PScalar3d<std::common_type_t<T, U>>(a.x * B.x + a.y * B.y + a.z * B.z);
}

// wedge product between a bivector and a vector (returns a trivector in 3d)
// wdg(A,b) = gr3( gpr(A,b) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& A, Vec3d<U> const& b)
{
    return PScalar3d<std::common_type_t<T, U>>(A.x * b.x + A.y * b.y + A.z * b.z);
}

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> and BiVec3d<T> projections and rejections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v1,
                                                              Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return dot(v1, v2) * Vec3d<ctype>(inv(v2));
}

// projection of v1 onto v2 (v2 must already be unitized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto_unitized(Vec3d<T> const& v1,
                                                                       Vec3d<U> const& v2)
{
    return dot(v1, v2) * v2; // v2 is already its own reverse
}

// projection of a vector v1 onto a bivector v2
// v_parallel = gpr( dot(v1,v2), inv(v2) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v1,
                                                              BiVec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    Vec3d<ctype> a = dot(v1, v2);
    BiVec3d<ctype> Bi = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // gpr(a,Bi) = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(gpr(a,Bi)) = dot(a,Bi)
    return Vec3d<ctype>(dot(a, Bi));
}

// projection of a vector v1 onto a unitized bivector v2
// u_parallel = gr1(gpr( dot(v1,v2), inv(v2) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
project_onto_unitized(Vec3d<T> const& v1, BiVec3d<U> const& v2)
{
    // requires v2 to be unitized

    using ctype = std::common_type_t<T, U>;
    Vec3d<ctype> a = dot(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec3d<ctype> Bi = -v2;
    // use the formular equivalent to the geometric product to save computational cost
    // gpr(a,Bi) = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(gpr(a,Bi)) = dot(a,Bi)
    return Vec3d<ctype>(dot(a, Bi));
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(gpr( wdg(v1,v2), inv(v2) ))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v1,
                                                             Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    BiVec3d<ctype> B = wdg(v1, v2);
    Vec3d<ctype> v2_inv = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // gpr(B,b_inv) = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(gpr(B,b_inv)) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    return Vec3d<ctype>(dot(B, v2_inv));
}

// rejection of vector v1 from a unitized vector v2
// v_perp = gr1(gpr( wdg(v1,v2), inv(v2) ))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from_unitized(Vec3d<T> const& v1,
                                                                      Vec3d<U> const& v2)
{
    // requires v2 to be unitized
    using ctype = std::common_type_t<T, U>;
    BiVec3d<ctype> B = wdg(v1, v2);
    Vec3d<ctype> v2_inv = v2; // v2 is its own inverse, if unitized
    // use the formular equivalent to the geometric product to save computational cost
    // gpr(B,b_inv) = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(gpr(B,b_inv)) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    return Vec3d<ctype>(dot(B, v2_inv));
}

// rejection of vector v1 from a bivector v2
// u_perp = gr1(gpr( wdg(v1,v2), inv(v2) ))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v1,
                                                             BiVec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    PScalar3d<ctype> a = wdg(v1, v2);
    BiVec3d<ctype> B = inv(v2);
    // trivector * bivector = vector
    return gpr(a, B);
}

// rejection of vector v1 from a unitized bivector v2
// u_perp = gpr( wdg(v1,v2), inv(v2) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
reject_from_unitized(Vec3d<T> const& v1, BiVec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    PScalar3d<ctype> a = wdg(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec3d<ctype> B = -v2;
    // trivector * bivector = vector (derived from full geometric product to save
    // costs)
    return Vec3d<ctype>(-a * B.x, -a * B.y, -a * B.z);
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// geometric product AB for fully populated 3d multivector
// gpr() ... geometric product
// Expensive! - Don't use if you don't have to! (64x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector vector = scalar + bivector)
// Ab = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector vector = vector + trivector)
// aB = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)  (vector bivector = vector + trivector)
// => see overloaded versions of gpr
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> gpr(MVec3d<T> const& A,
                                                      MVec3d<U> const& B)
{
    // geometric product with a fully populated 3d multivector
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 -
               A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 - A.c4 * B.c7 -
               A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c4;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 + A.c4 * B.c3 -
               A.c5 * B.c7 - A.c6 * B.c1 - A.c7 * B.c5;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
               A.c5 * B.c1 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0 -
               A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c6 +
               A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c7 - A.c4 * B.c5 +
               A.c5 * B.c4 + A.c6 * B.c0 + A.c7 * B.c3;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
               A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product ab between two vectors (returns a multivector of the even subalgebra)
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector vector = scalar + bivector)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> gpr(Vec3d<T> const& a,
                                                        Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(dot(a, b)), wdg(a, b));
}

// define geometric multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product Ab for a bivector and a vector (returns a multivector)
// Ab = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)
// => bivector vector = vector + trivector
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> gpr(BiVec3d<T> const& A,
                                                        Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(dot(A, b), wdg(A, b));
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& A,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, b);
}

// geometric product aB for a vector and a bivector (returns a multivector)
// aB = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)
// => bivector vector = vector + trivector
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> gpr(Vec3d<T> const& a,
                                                        BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(dot(a, B), wdg(a, B));
}

// define geometric multiplication with operator*(a,B) as an alias for gpr(a,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, B);
}

// geometric product ab between two bivectors (returns an even grade multivector)
// ab = gr0(ab) + gr2(ab)
// => bivector bivector = scalar + bivector, in 3D)
//
// the full geometric bivector product only exists in 4 dimensional spaces:
// ab = a*b + axb + a^b = gr0(ab) + gr2(ab) + gr4(ab)
// In 3D we don't have gr4(ab) and thus only the terms up to grade 3 remain.
// The bivector product axb = 0.5*(ab-ba)is called the commutator product.
//
// ab = dot(a,b) + cmt(a,b) + wgd(a,b)  (in 4D and higher dimensional spaces)
// ab = dot(a,b) + cmt(a,b)             (in 3D)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> gpr(BiVec3d<T> const& a,
                                                        BiVec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(dot(a, b)), cmt(a, b));
}

// define geometric multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& a,
                                                              BiVec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product AB of a trivector A multiplied from the left
// to the multivector B
// gpr(trivector, multivector) => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> gpr(PScalar3d<T> A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MVec3d<ctype>(-B.c7, -B.c4, -B.c5, -B.c6, B.c1, B.c2, B.c3, B.c0);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                            MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product Ab of a trivector A multiplied from the left
// to the vector b
// gpr(trivector, vector) => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> gpr(PScalar3d<T> A, Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * BiVec3d<ctype>(b.x, b.y, b.z);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                             Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, b);
}

// geometric product AB of a trivector A multiplied from the left
// to the multivector from the uneven subalgebra B
// gpr(trivector, uneven multivector) => returns a even multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> gpr(PScalar3d<T> A,
                                                        MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) *
           MVec3d_E<ctype>(Scalar3d<ctype>(-B.c3), BiVec3d<ctype>(B.c0, B.c1, B.c2));
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product Ab of a trivector A multiplied from the left
// to the bivector b
// gpr(trivector, bivector) => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> gpr(PScalar3d<T> A, BiVec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * Vec3d<ctype>(b.x, b.y, b.z);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                           BiVec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, b);
}

// geometric product AB of a trivector A multiplied from the left
// to the multivector from the even subalgebra B
// gpr(trivector, even multivector) => returns an uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> gpr(PScalar3d<T> A,
                                                        MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) *
           MVec3d_U<ctype>(Vec3d<ctype>(-B.c1, -B.c2, -B.c3), PScalar3d<ctype>(B.c0));
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product AB of a trivector B multiplied from the right
// to the multivector A
// gpr(multivector, trivector) => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> gpr(MVec3d<T> const& A, PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(-A.c7, -A.c4, -A.c5, -A.c6, A.c1, A.c2, A.c3, A.c0) * ctype(B);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product aB of a trivector B multiplied from the right
// to the vector a
// gpr(vector, trivector) => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> gpr(Vec3d<T> const& a, PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(a.x, a.y, a.z) * ctype(B);
}

// define geometric multiplication with operator*(a,B) as an alias for gpr(a,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                             PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, B);
}

// geometric product AB of an uneven multivector A multiplied from the right
// by the trivector B (=pseudoscalar in 3d)
// gpr(uneven multivector, trivector) => returns a even multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> gpr(MVec3d_U<U> const& A,
                                                        PScalar3d<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-A.c3), BiVec3d<ctype>(A.c0, A.c1, A.c2)) *
           ctype(B);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<U> const& A,
                                                              PScalar3d<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product aB of a trivector B multiplied from the right
// to the bivector a
// gpr(bivector, trivector) => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> gpr(BiVec3d<T> const& a, PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return -Vec3d<ctype>(a.x, a.y, a.z) * ctype(B);
}

// define geometric multiplication with operator*(a,B) as an alias for gpr(a,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& a,
                                                           PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, B);
}

// geometric product AB of an even multivector A multiplied from the right
// by the trivector B (=pseudoscalar in 3d)
// gpr(even multivector, trivector) => returns an uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> gpr(MVec3d_E<U> const& A,
                                                        PScalar3d<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-A.c1, -A.c2, -A.c3), PScalar3d<ctype>(A.c0)) *
           ctype(B);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<U> const& A,
                                                              PScalar3d<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// geometric product AB of two trivectors
// gpr(trivector, trivector) => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> gpr(PScalar3d<T> A, PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return -ctype(A) * ctype(B); // trivectors square to -1
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> operator*(PScalar3d<T> A, PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

// return conjugate complex of a trivector
// i.e. the reverse in nomenclature of multivectors
template <typename T> inline PScalar3d<T> rev(PScalar3d<T> A)
{
    // the 3d trivector switches sign on reversion
    return PScalar3d<T>(-T(A));
}

////////////////////////////////////////////////////////////////////////////////
// MVec3d_E<T> / MVec3d_U<T> operations
////////////////////////////////////////////////////////////////////////////////

// geometric product ab for two multivectors from the even subalgebra (3d case)
// a  = gr0(a)  + gr2(a)  (scalar + bivector)
// b  = gr0(b)  + gr2(b)  (scalar + bivector)
// ab = gr0(ab) + gr2(ab) (scalar + bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> gpr(MVec3d_E<T> const& a,
                                                        MVec3d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(
        Scalar3d<ctype>(a.c0 * b.c0 - a.c1 * b.c1 - a.c2 * b.c2 - a.c3 * b.c3),
        BiVec3d<ctype>(a.c0 * b.c1 + a.c1 * b.c0 - a.c2 * b.c3 + a.c3 * b.c2,
                       a.c0 * b.c2 + a.c1 * b.c3 + a.c2 * b.c0 - a.c3 * b.c1,
                       a.c0 * b.c3 - a.c1 * b.c2 + a.c2 * b.c1 + a.c3 * b.c0));
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& a,
                                                              MVec3d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product ab for two multivectors from the uneven subalgebra (3d case)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> gpr(MVec3d_U<T> const& a,
                                                        MVec3d_U<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(
        Scalar3d<ctype>(a.c0 * b.c0 + a.c1 * b.c1 + a.c2 * b.c2 - a.c3 * b.c3),
        BiVec3d<ctype>(a.c0 * b.c3 + a.c1 * b.c2 - a.c2 * b.c1 + a.c3 * b.c0,
                       -a.c0 * b.c2 + a.c1 * b.c3 + a.c2 * b.c0 + a.c3 * b.c1,
                       a.c0 * b.c1 - a.c1 * b.c0 + a.c2 * b.c3 + a.c3 * b.c2));
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& a,
                                                              MVec3d_U<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product Ab of a multivector A from the even subalgebra (3d case)
// with a vector b from the right
// => returns an uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> gpr(MVec3d_E<T> const& A,
                                                        Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(A.c0 * b.x - A.c2 * b.z + A.c3 * b.y,
                                        A.c0 * b.y + A.c1 * b.z - A.c3 * b.x,
                                        A.c0 * b.z - A.c1 * b.y + A.c2 * b.x),
                           PScalar3d<ctype>(A.c1 * b.x + A.c2 * b.y + A.c3 * b.z));
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, b);
}

// geometric product ab of a multivector a from the even subalgebra (3d case)
// with a bivector b from the right
// => returns an even multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> gpr(MVec3d_E<T> const& a,
                                                        BiVec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-a.c1 * b.x - a.c2 * b.y - a.c3 * b.z),
                           BiVec3d<ctype>(a.c0 * b.x - a.c2 * b.z + a.c3 * b.y,
                                          a.c0 * b.y + a.c1 * b.z - a.c3 * b.x,
                                          a.c0 * b.z - a.c1 * b.y + a.c2 * b.x));
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& a,
                                                              BiVec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product aB a multivector B from the even subalgebra (3d case)
// with a vector from the left
// => returns an uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> gpr(Vec3d<T> const& a,
                                                        MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(a.x * B.c0 - a.y * B.c3 + a.z * B.c2,
                                        a.x * B.c3 + a.y * B.c0 - a.z * B.c1,
                                        -a.x * B.c2 + a.y * B.c1 + a.z * B.c0),
                           PScalar3d<ctype>(a.x * B.c1 + a.y * B.c2 + a.z * B.c3));
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, B);
}

// geometric product ab of a multivector b from the even subalgebra (3d case)
// with a bivector a from the right
// => returns an even multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> gpr(BiVec3d<T> const& a,
                                                        MVec3d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-a.x * b.c1 - a.y * b.c2 - a.z * b.c3),
                           BiVec3d<ctype>(a.x * b.c0 - a.y * b.c3 + a.z * b.c2,
                                          a.x * b.c3 + a.y * b.c0 - a.z * b.c1,
                                          -a.x * b.c2 + a.y * b.c1 + a.z * b.c0));
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& a,
                                                              MVec3d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

// geometric product AB of a multivector A from the uneven subalgebra (3d case)
// with a multivector B of the even subalgebra
// => returns an uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> gpr(MVec3d_U<T> const& A,
                                                        MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(
        Vec3d<ctype>(A.c0 * B.c0 - A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1,
                     A.c0 * B.c3 + A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2,
                     -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar3d<ctype>(A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0));
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                              MVec3d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, b);
}

// geometric product AB of a multivector A from the even subalgebra (3d case)
// with a multivector B of the uneven subalgebra
// => returns an uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> gpr(MVec3d_E<T> const& A,
                                                        MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(
        Vec3d<ctype>(A.c0 * B.c0 - A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1,
                     A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c3 - A.c3 * B.c0,
                     A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar3d<ctype>(A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2));
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(A, B);
}

//////////////////////////////////////////////////////////////////////////////////////////
// exponential function with bivector as argument for setup of quaternions
// as geometric multivector with a scalar and a bivector part
// MVec3d_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
//
// quaternion: q = a + b I with I being the bivector in brackets above
//             representing a plane in the algebra G^3
//
//             a rotation in 3D is represented by the plane and the
//             size of the rotation, the later is given by the angle
//             theta, which is the magnitude of the bivector
//
// Inputs:
//         - an arbitray bivector representing the oriented plane of rotation
//           (does not need to be unitized)
//         - a rotation angle
// Output:
//         - a rotor representing the rotation
//
// HINT:     For a rotation around an axis n (n = unitized(Vec3d<T>))
//           use the bivector B = n*I_3d  => B = gpr(Vec3d<T>,PScalar3d<T>
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> exp(BiVec3d<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(std::cos(theta)),
                           unitized(I) * std::sin(theta));
}

//////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//       1.) an arbitray bivector representing the oriented plane of rotation
//           (does not need to be unitized, defines what is a posive rotation angle)
//       2.) a rotation angle in that plane
// Output:
//           a rotor representing the requested rotation,
//           for applying the sandwich product as in rotate(v,rotor)
//
//////////////////////////////////////////////////////////////////////////////////////////
//
// for a rotation about an axis n (n = unitized vector) choose the ansatz n*B = I_3d
// an multiply both sides with n from the left (remember n*n = |n|^2 = 1)
//
// => choose: B = n*I_3d
//
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> rotor(BiVec3d<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec3d_E<ctype>(Scalar3d<ctype>(std::cos(half_angle)),
                           unitized(I) * std::sin(half_angle));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rotate(Vec3d<T> const& v,
                                                        MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> reverse_rotor = rev(rotor);
    // MVec3d_U<ctype> tmp = rotor * v;
    // MVec3d_U<ctype> res = tmp * reverse rotor;

    // trivector part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    return Vec3d<ctype>(gr1<ctype>(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> rotate(BiVec3d<T> const& v,
                                                          MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> reverse_rotor = rev(rotor);
    // MVec3d_E<ctype> tmp = rotor * v;
    // MVec3d_E<ctype> res = tmp * reverse rotor;

    // scalar part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    return BiVec3d<ctype>(gr2<ctype>(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> rotate(MVec3d<T> const& v,
                                                         MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(rotor * v * rev(rotor));
}

// return the dual(M) of the multivector M
// if M represents the subspace B as subspace of R^2 then
// dual(M) represents the orthogonal subspace B^perp (perpendicular to B)
// => returns the orthogonal complement
//
// dual by left multiplication with Im_3d
// as defined in Doran/Lasenby "GA for physicists"
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> dual(MVec3d<T> const& M)
{
    return MVec3d<T>(-M.c7, -M.c4, -M.c5, -M.c6, M.c1, M.c2, M.c3, M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> dual(MVec3d_E<T> const& M)
{
    return MVec3d_U<T>(-M.c1, -M.c2, -M.c3, M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> dual(MVec3d_U<T> const& M)
{
    return MVec3d_E<T>(-M.c3, M.c0, M.c1, M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> dual(PScalar3d<T> const& ps)
{
    return Scalar3d<T>(-T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> dual(Vec3d<T> const& v)
{
    return BiVec3d<T>(v.x, v.y, v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> dual(BiVec3d<T> const& B)
{
    return Vec3d<T>(-B.x, -B.y, -B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual(Scalar3d<T> const& s)
{
    return PScalar3d<T>(T(s));
}

} // namespace hd::ga
