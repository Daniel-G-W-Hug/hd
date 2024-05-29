#pragma once

// author: Daniel Hug, 2024

#include <algorithm> // std::clamp
#include <cmath>     // std::abs, std::sin, std::cos
#include <concepts>  // std::floating_point<T>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "hd_ga_cfg_value_t.hpp"

#include "hd_ga_cfg_vec2d.hpp"

#include "hd_ga_cfg_mvec2d.hpp"
#include "hd_ga_cfg_mvec2d_e.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of v1 onto v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto(Vec2d<T> const& v1,
                                                              Vec2d<U> const& v2)
{
    return dot(v1, v2) * inv(v2);
}

// projection of v1 onto v2 (v2 must already be unitized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> project_onto_unitized(Vec2d<T> const& v1,
                                                                       Vec2d<U> const& v2)
{
    // requires v2 to be unitized
    return dot(v1, v2) * v2;
}

// rejection of v1 from v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from(Vec2d<T> const& v1,
                                                             Vec2d<U> const& v2)
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
inline constexpr Vec2d<std::common_type_t<T, U>> reject_from_unitized(Vec2d<T> const& v1,
                                                                      Vec2d<U> const& v2)
{
    // requires v2 to be unitized

    using ctype = std::common_type_t<T, U>;
    // version using geometric algebra wedge product manually computed
    // from "wdg(v1,v2)*inv(v2)" + v2 being already it's own inverse
    PScalar2d<ctype> w = wdg(v1, v2); // bivector with component e12
    return Vec2d<ctype>(v2.y * w, -v2.x * w);
}

// reflect a vector u on a hyperplane B orthogonal to vector b
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n
// (e.g. a line through the origin in 2d space)
// orthogonal to vector b: the hyperplane is dual to b
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reflect_on_hyp(Vec2d<T> const& u,
                                                                Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-b * u * inv(b));
}

// reflect a vector u another vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> reflect_on_vec(Vec2d<T> const& u,
                                                                Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(b * u * inv(b));
}

////////////////////////////////////////////////////////////////////////////////
// MVec2d<T> geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product ab for fully populated 2d multivector
// gpr() ... geometric product
// Expensive! - Don't use if you don't have to! (16x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors, e.g.:
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab) (vector vector = scalar + bivector)
//
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
    return gpr(A, B);
}

// geometric product ab for two vectors (returns a multivector of the even subalgebra)
// ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
// => vector vector = scalar + bivector
//
// HINT: if a full 2d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> gpr(Vec2d<T> const& a,
                                                        Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(Scalar<ctype>(dot(a, b)), wdg(a, b));
}

// define geometric multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                              Vec2d<U> const& b)
{
    return gpr(a, b);
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
    return gpr(A, B);
}

// geometric product AB of a bivector A multiplied from the left
// to the multivector from the even subalgebra B (MVec2d_E)
// gpr(bivector, even grade multivector) => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> gpr(PScalar2d<T> A,
                                                        MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MVec2d_E<ctype>(-B.c1, B.c0);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(PScalar2d<T> A,
                                                              MVec2d_E<U> const& B)
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
    return gpr(A, b);
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
    return gpr(A, B);
}

// geometric product AB of a bivector B multiplied from the right
// to the even grade multivector A (MVec2d_E)
// gpr(even grade multivector, bivector) => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> gpr(MVec2d_E<T> const& A,
                                                        PScalar2d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(-A.c1, A.c0) * ctype(B);
}

// define geometric multiplication with operator*(A,B) as an alias for gpr(A,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
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
    return gpr(a, B);
}

// geometric product AB of two bivectors (pseudoscalars in 2d)
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
    return gpr(A, B);
}

// geometric product aB for a full 2d multivector B with a vector a
// => return a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(Vec2d<T> const& a,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(a.x * B.c1 + a.y * B.c2, a.x * B.c0 - a.y * B.c3,
                         a.x * B.c3 + a.y * B.c0, a.x * B.c2 - a.y * B.c1);
}

// define geometric multiplication with operator*(a,B) as an alias for gpr(a,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                            MVec2d<U> const& B)
{
    return gpr(a, B);
}

// geometric product aB for a full 2d multivector B with a multivector from
// the even subalgebra a (only gr0 and gr2 components, i.e. a MVec2d_E)
// => product with MVec2d_E from the left
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MVec2d_E<T> const& a,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(a.c0 * B.c0 - a.c1 * B.c3, a.c0 * B.c1 + a.c1 * B.c2,
                         a.c0 * B.c2 - a.c1 * B.c1, a.c0 * B.c3 + a.c1 * B.c0);
}

// define geometric multiplication with operator*(a,B) as an alias for gpr(a,B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& a,
                                                            MVec2d<U> const& B)
{
    return gpr(a, B);
}

// geometric product ab for a 2d vector b with a multivector from
// the even subalgebra a (only gr0 and gr2 components, i.e. a MVec2d_E)
// => product with MVec2d_E from the left
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> gpr(MVec2d_E<T> const& a,
                                                     Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(a.c0 * b.x + a.c1 * b.y, a.c0 * b.y - a.c1 * b.x);
}

// define geometric multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& a,
                                                           Vec2d<U> const& b)
{
    return gpr(a, b);
}

// geometric product Ab for a full 2d multivector A with a multivector from
// the even subalgebra b (only gr0 and gr2 components, i.e. a MVec2d_E)
// => product with MVec2d_E from the right
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MVec2d<T> const& A,
                                                      MVec2d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c0 * b.c0 - A.c3 * b.c1, A.c1 * b.c0 - A.c2 * b.c1,
                         A.c1 * b.c1 + A.c2 * b.c0, A.c0 * b.c1 + A.c3 * b.c0);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            MVec2d_E<U> const& b)
{
    return gpr(A, b);
}

// geometric product Ab for a full 2d multivector A with a vector b
// => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> gpr(MVec2d<T> const& A,
                                                      Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c1 * b.x + A.c2 * b.y, A.c0 * b.x + A.c3 * b.y,
                         -A.c3 * b.x + A.c0 * b.y, -A.c2 * b.x + A.c1 * b.y);
}

// define geometric multiplication with operator*(A,b) as an alias for gpr(A,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                            Vec2d<U> const& b)
{
    return gpr(A, b);
}

// geometric product aB for a full 2d vector a with a multivector from
// the even subalgebra b (only gr0 and gr2 components, i.e. a MVec2d_E)
// => product with MVec2d_E from the right
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> gpr(Vec2d<T> const& a,
                                                     MVec2d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(a.x * b.c0 - a.y * b.c1, a.x * b.c1 + a.y * b.c0);
}

// define geometric multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& a,
                                                           MVec2d_E<U> const& b)
{
    return gpr(a, b);
}

// geometric product ab for two multivectors from the even subalgebra (2d case)
// a  = gr0(a)  + gr2(a)  (scalar + bivector)
// b  = gr0(b)  + gr2(b)  (scalar + bivector)
// ab = gr0(ab) + gr2(ab) (scalar + bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> gpr(MVec2d_E<T> const& a,
                                                        MVec2d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(a.c0 * b.c0 - a.c1 * b.c1, a.c0 * b.c1 + a.c1 * b.c0);
}

// define complex multiplication with operator*(a,b) as an alias for gpr(a,b)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& a,
                                                              MVec2d_E<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return gpr<ctype>(a, b);
}

////////////////////////////////////////////////////////////////////////////////
// 2d rotation operations
////////////////////////////////////////////////////////////////////////////////

// exponential function for setup of complex numbers and rotations
// as geometric multivectors with a scalar and a bivector part
//
// r = 1 is the vector length of the complex number in polar form
// theta is the bivector angle (i.e. a multiple of the bivector I_2d)
// such that uv = r exp(I_2d, theta) = a + I_2d b
// with r = |u| |v| = sqrt(a^2 + b^2) = 1
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> exp([[maybe_unused]] PScalar2d<T> I,
                                                        U theta)
{
    // PScalar2d<T> is just used here for a unique overload of exp()
    // and to make the function signature similar to the 3D case
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(Scalar<ctype>(std::cos(theta)),
                           PScalar2d<ctype>(std::sin(theta)));
}

//////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//         - a 2d pseudoscalar representing the plane of 2d space
//         - a rotation angle in the plane of 2d space
// Output:
//         - a rotor representing the requested rotation,
//           when applying the sandwich product with the rotor as in rotate(v,rotor)
//
//////////////////////////////////////////////////////////////////////////////////////////
//
// implemented here to make it formally the same with the 3d case (and potentially higher
// dimensional applications). In 2d the rotation can be directly expressed with less
// effort as
//
// exp(I_2d, -theta) * v = v * exp(I_2d, theta)  to express a 2d rotation of vector v by
//                                               the angle theta
//
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> rotor([[maybe_unused]] PScalar2d<T> I,
                                                          U theta)
{
    // PScalar2d<T> is just used here to be able to overload exp
    // and to make the function similar to the 3D case
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec2d_E<ctype>(Scalar<ctype>(std::cos(half_angle)),
                           PScalar2d<ctype>(std::sin(half_angle)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec2d<std::common_type_t<T, U>> rotate(Vec2d<T> const& v,
                                                        MVec2d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(rotor * v * rev(rotor));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d<std::common_type_t<T, U>> rotate(MVec2d<T> const& M,
                                                         MVec2d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(rotor * M * rev(rotor));
}

////////////////////////////////////////////////////////////////////////////////
// 2d duality operations
////////////////////////////////////////////////////////////////////////////////

// return the dual(M) of the multivector M
// if M represents the subspace B as subspace of R^2 then
// dual(M) represents the orthogonal subspace B^perp (perpendicular to B)
// => returns the orthogonal complement
//
// dual by left multiplication with Im_2d
// as defined in Doran/Lasenby "GA for physicists"
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar<T> dual2d(PScalar2d<T> ps)
{
    return Scalar<T>(-T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> dual2d(Scalar<T> s)
{
    return PScalar2d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> dual2d(Vec2d<T> const& v)
{
    return Vec2d<T>(v.y, -v.x);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> dual2d(MVec2d_E<T> const& M)
{
    return MVec2d_E<T>(-M.c1, M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> dual2d(MVec2d<T> const& M)
{
    return MVec2d<T>(-M.c3, M.c2, -M.c1, M.c0);
}

} // namespace hd::ga
