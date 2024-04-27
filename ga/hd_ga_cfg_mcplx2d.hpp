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
// MCplx2d<T> M = c0 + c1 I (with I being the bivector of the plane e1^e2)
//
// This is the definition of a multivector in the even subalgebra of G<2,0,0>,
// i.e. it models only multivectors with even grades 0 and 2 in the plane e1^e2
// (complex numbers).
// This subalgebra is closed under addition and multiplication.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.


template <typename T = value_t>
    requires(std::floating_point<T>)
struct MCplx2d {

    // ctors

    // (all grades = 0)
    MCplx2d() = default;

    // assign all components
    MCplx2d(T s, T ps) : c0(s), c1(ps) {}

    // assign a scalar part exclusively (other grades = 0)
    MCplx2d(Scalar<T> s) : c0(s) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MCplx2d(PScalar2d<T> ps) : c1(ps) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to full geometric product)
    MCplx2d(Scalar<T> s, PScalar2d<T> ps) : c0(s), c1(ps) {}

    // floating point type conversion
    template <typename U>
        requires(std::floating_point<U>)
    MCplx2d(MCplx2d<U> const& v) : c0(v.c0), c1(v.c1)
    {
    }


    T c0{}; // scalar component
    T c1{}; // bivector component (2d pseudoscalar)

    // equality
    template <typename U>
        requires(std::floating_point<U>)
    bool operator==(MCplx2d<U> const& rhs) const
    {
        // componentwise comparison
        // equality implies same magnitude and direction
        // comparison is not exact, but accepts epsilon deviations
        auto abs_delta_c0 = std::abs(rhs.c0 - c0);
        auto abs_delta_c1 = std::abs(rhs.c1 - c1);
        auto constexpr delta_eps =
            std::common_type_t<T, U>(5.0) *
            std::max<std::common_type_t<T, U>>(std::numeric_limits<T>::epsilon(),
                                               std::numeric_limits<U>::epsilon());
        if (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps) return true;
        return false;
    }

    // unary minus (must be declared a friend otherwise doesn't work)
    friend inline constexpr MCplx2d<T> operator-(MCplx2d<T> const& v)
    {
        return MCplx2d<T>(-v.c0, -v.c1);
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, MCplx2d<U> const& v);
};

// for printing via iostream
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, MCplx2d<T> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////
// MCplx2d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// add multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator+(MCplx2d<T> const& v1,
                                                             MCplx2d<U> const& v2)
{
    return MCplx2d<std::common_type_t<T, U>>(v1.c0 + v2.c0, v1.c1 + v2.c1);
}

// substract multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator-(MCplx2d<T> const& v1,
                                                             MCplx2d<U> const& v2)
{
    return MCplx2d<std::common_type_t<T, U>>(v1.c0 - v2.c0, v1.c1 - v2.c1);
}


// multiply a multivector with a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator*(MCplx2d<T> const& v, U s)
{
    return MCplx2d<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator*(T s, MCplx2d<U> const& v)
{
    return MCplx2d<std::common_type_t<T, U>>(v.c0 * s, v.c1 * s);
}

// devide a multivector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MCplx2d<std::common_type_t<T, U>> operator/(MCplx2d<T> const& v, U s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return MCplx2d<std::common_type_t<T, U>>(v.c0 * inv, v.c1 * inv);
}

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector (= pseudoscalar in 2d)

template <typename T> inline constexpr Scalar<T> gr0(MCplx2d<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr PScalar2d<T> gr2(MCplx2d<T> const& v)
{
    return PScalar2d<T>(v.c1);
}

} // namespace hd::ga


// ////////////////////////////////////////////////////////////////////////////////
// // printing support via fmt library
// ////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MCplx2d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MCplx2d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{})", nested(v.c0), nested(v.c1));
    }
};
