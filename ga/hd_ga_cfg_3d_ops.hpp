#pragma once

// author: Daniel Hug, 2024

#include <cmath>    // abs, sqrt, acos
#include <concepts> // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "hd_ga_cfg_bivec3d.hpp"
#include "hd_ga_cfg_value_t.hpp"
#include "hd_ga_cfg_vec3d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> & BiVec3d<T> geometric operations
////////////////////////////////////////////////////////////////////////////////

// return the dot product of two vectors (= a scalar)
// coordinate free definition: dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(const Vec3d<T>& v1, const Vec3d<U>& v2)
{
    // this implementation is only valid in an orthonormal basis
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// return squared magnitude of vector
template <typename T> inline T sq_nrm(const Vec3d<T>& v) { return dot(v, v); }

// return magnitude of vector
template <typename T> inline T nrm(const Vec3d<T>& v) { return std::sqrt(dot(v, v)); }

// return the dot product of a bivector and a vector (= a vector)
// dot(A,b) = gr0( gpr(A,b) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> dot(const BiVec3d<T>& v1, const Vec3d<U>& v2)
{
    // this implementation is only valid in an orthonormal basis
    return Vec3d<std::common_type_t<T, U>>(
        v1.z * v2.y - v1.y * v2.z, v1.x * v2.z - v1.z * v2.x, v1.y * v2.x - v1.x * v2.y);
}

// return the dot product of a vector and a bivector (= a vector)
// dot(a,B) = gr1( gpr(a,B) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline Vec3d<std::common_type_t<T, U>> dot(const Vec3d<T>& v1, const BiVec3d<U>& v2)
{
    // this implementation is only valid in an orthonormal basis
    return Vec3d<std::common_type_t<T, U>>(
        v1.z * v2.y - v1.y * v2.z, v1.x * v2.z - v1.z * v2.x, v1.y * v2.x - v1.x * v2.y);
}

// return dot product of two bivectors A and B (= a scalar)
// dot(A,B) = gr0( gpr(A,B) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> dot(const BiVec3d<T>& v1, const BiVec3d<U>& v2)
{
    // this implementation is only valid in an orthonormal basis
    return -v1.x * v2.x - v1.y * v2.y - v1.z * v2.z;
}

// return squared magnitude of bivector
//
// TODO: Check whether this the right way to calculate the magnitude
//
template <typename T> inline T sq_nrm(const BiVec3d<T>& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

// return magnitude of bivector
template <typename T> inline T nrm(const BiVec3d<T>& v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// return a vector unitized to nrm(v) == 1.0
template <typename T> inline Vec3d<T> unitized(const Vec3d<T>& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return Vec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return a bivector unitized to nrm(v) == 1.0
template <typename T> inline BiVec3d<T> unitized(const BiVec3d<T>& v)
{
    T n = nrm(v);
    if (n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for normalization" +
                                 std::to_string(n) + "\n");
    }
    T inv = 1.0 / n; // for multiplication with inverse of norm
    return BiVec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the multiplicative inverse of the vector
template <typename T> inline Vec3d<T> inv(const Vec3d<T>& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = 1.0 / sq_n; // inverse of squared norm for a vector
    return Vec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the multiplicative inverse of the bivector
template <typename T> inline BiVec3d<T> inv(const BiVec3d<T>& v)
{
    T sq_n = sq_nrm(v);
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion" +
                                 std::to_string(sq_n) + "\n");
    }
    T inv = -1.0 / sq_n; // negative inverse of squared norm for a bivector
    return BiVec3d<T>(v.x * inv, v.y * inv, v.z * inv);
}

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(const Vec3d<T>& v1, const Vec3d<U>& v2)
{
    std::common_type_t<T, U> nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<std::common_type_t<T, U>>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    return std::acos(dot(v1, v2) / nrm_prod);
}

// return the angle between of a vector and a bivector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(const Vec3d<T>& v1, const BiVec3d<U>& v2)
{
    std::common_type_t<T, U> nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<std::common_type_t<T, U>>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    return std::acos(std::abs(dot(v1, v2)) / nrm_prod);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(const BiVec3d<T>& v1, const Vec3d<U>& v2)
{
    std::common_type_t<T, U> nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<std::common_type_t<T, U>>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    return std::acos(std::abs(dot(v1, v2)) / nrm_prod);
}

// return the angle between two bivectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(const BiVec3d<T>& v1, const BiVec3d<U>& v2)
{
    std::common_type_t<T, U> nrm_prod = nrm(v1) * nrm(v2);
    if (nrm_prod < std::numeric_limits<std::common_type_t<T, U>>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle" +
            std::to_string(nrm_prod) + "\n");
    }
    return std::acos(std::abs(dot(v1, v2)) / nrm_prod);
}

// cross-product between two vectors (returns a vector in 3d)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec3d<std::common_type_t<T, U>> cross(const Vec3d<T>& v1, const Vec3d<U>& v2)
{
    return Vec3d<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// wedge product between two vectors (returns a bivector in 3d)
// coordinate-free definition: wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline BiVec3d<std::common_type_t<T, U>> wdg(const Vec3d<T>& v1, const Vec3d<U>& v2)
{
    return BiVec3d<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// wedge product between a vector and a bivector (returns a trivector in 3d)
// wdg(a,B) = gr3( gpr(a,B) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3d_t<std::common_type_t<T, U>> wdg(const Vec3d<T>& v1, const BiVec3d<U>& v2)
{
    return PScalar3d_t<std::common_type_t<T, U>>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// wedge product between a bivector and a vector (returns a trivector in 3d)
// wdg(A,b) = gr3( gpr(A,b) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline PScalar3d_t<std::common_type_t<T, U>> wdg(const BiVec3d<T>& v1, const Vec3d<U>& v2)
{
    return PScalar3d_t<std::common_type_t<T, U>>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// projection of a vector v1 onto vector v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(const Vec3d<T>& v1,
                                                              const Vec3d<U>& v2)
{
    return dot(v1, v2) * Vec3d<std::common_type_t<T, U>>(inv(v2));
}

// projection of v1 onto v2 (v2 must already be unitized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto_unitized(const Vec3d<T>& v1,
                                                                       const Vec3d<U>& v2)
{
    return dot(v1, v2) * v2;
}

// projection of a vector v1 onto a bivector v2
// v_parallel = gpr( dot(v1,v2), inv(v2) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(const Vec3d<T>& v1,
                                                              const BiVec3d<U>& v2)
{
    using ctype = std::common_type_t<T, U>;
    Vec3d<ctype> a = dot(v1, v2);
    BiVec3d<ctype> Bi = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // aBi = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(aBi) = dot(a,Bi)
    return Vec3d<std::common_type_t<T, U>>(dot(a, Bi));
}

// projection of a vector v1 onto a unitized bivector v2
// u_parallel = gr1(gpr( dot(v1,v2), inv(v2) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
project_onto_unitized(const Vec3d<T>& v1, const BiVec3d<U>& v2)
{
    // requires v2 to be unitized

    Vec3d<std::common_type_t<T, U>> a = dot(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec3d<std::common_type_t<T, U>> Bi = -v2;
    // use the formular equivalent to the geometric product to save computational cost
    // aBi = dot(a,Bi) + wdg(a,Bi)
    // v_parallel = gr1(aBi) = dot(a,Bi)
    return Vec3d<std::common_type_t<T, U>>(dot(a, Bi));
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(gpr( wdg(v1,v2), inv(v2) ))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(const Vec3d<T>& v1,
                                                             const Vec3d<U>& v2)
{
    BiVec3d<std::common_type_t<T, U>> A = wdg(v1, v2);
    Vec3d<std::common_type_t<T, U>> bi = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // Abi = dot(A,bi) + wdg(A,bi)
    // v_perp = gr1(Abi) = dot(A,bi)
    return Vec3d<std::common_type_t<T, U>>(dot(A, bi));
}

// rejection of vector v1 from a unitized vector v2
// v_perp = gr1(gpr( wdg(v1,v2), inv(v2) ))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from_unitized(const Vec3d<T>& v1,
                                                                      const Vec3d<U>& v2)
{
    // requires v2 to be unitized

    BiVec3d<std::common_type_t<T, U>> A = wdg(v1, v2);
    Vec3d<std::common_type_t<T, U>> bi = v2;
    // use the formular equivalent to the geometric product to save computational cost
    // Abi = dot(A,bi) + wdg(A,bi)
    // v_perp = gr1(Abi) = dot(A,bi)
    return Vec3d<std::common_type_t<T, U>>(dot(A, bi));
}

// rejection of vector v1 from a bivector v2
// u_perp = gr1(gpr( wdg(v1,v2), inv(v2) ))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(const Vec3d<T>& v1,
                                                             const BiVec3d<U>& v2)
{
    PScalar3d_t<std::common_type_t<T, U>> a = wdg(v1, v2);
    BiVec3d<std::common_type_t<T, U>> B = inv(v2);
    // trivector * bivector = vector (derived from full geometric product to save
    // costs)
    return Vec3d<std::common_type_t<T, U>>(-a * B.x, -a * B.y, -a * B.z);
}

// rejection of vector v1 from a unitized bivector v2
// u_perp = gpr( wdg(v1,v2), inv(v2) )
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
reject_from_unitized(const Vec3d<T>& v1, const BiVec3d<U>& v2)
{
    PScalar3d_t<std::common_type_t<T, U>> a = wdg(v1, v2);
    // up to the sign v2 already is it's own inverse
    BiVec3d<std::common_type_t<T, U>> B = -v2;
    // trivector * bivector = vector (derived from full geometric product to save
    // costs)
    return Vec3d<std::common_type_t<T, U>>(-a * B.x, -a * B.y, -a * B.z);
}

} // namespace hd::ga
