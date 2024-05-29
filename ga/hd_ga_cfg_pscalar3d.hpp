#pragma once

// author: Daniel Hug, 2024

#include <concepts> // std::floating_point

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "hd_ga_cfg_value_t.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// PScalar3d<T> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> operator-(PScalar3d<T> v)
{
    return PScalar3d<T>(-value_t(v));
}

// adding pseudoscalars
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator+(PScalar3d<T> v1,
                                                               PScalar3d<U> v2)
{
    return PScalar3d<std::common_type_t<T, U>>(value_t(v1) + value_t(v2));
}

// substracting pseudoscalars
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator-(PScalar3d<T> v1,
                                                               PScalar3d<U> v2)
{
    return PScalar3d<std::common_type_t<T, U>>(value_t(v1) - value_t(v2));
}


// multiply a vector with a scalar (in both constellations)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator*(PScalar3d<T> v, U s)
{
    return PScalar3d<std::common_type_t<T, U>>(value_t(v) * s);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator*(T s, PScalar3d<U> v)
{
    return PScalar3d<std::common_type_t<T, U>>(s * value_t(v));
}

// devide a vector by a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator/(PScalar3d<T> v, U s)
{
    if (s == 0.0) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    U inv = 1.0 / s; // for multiplicaton with inverse value
    return PScalar3d<std::common_type_t<T, U>>(value_t(v) * inv);
}

////////////////////////////////////////////////////////////////////////////////
// PScalar3d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return squared magnitude of the pseudoscalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr T sq_nrm(PScalar3d<T> ps)
{
    return T(ps) * T(ps);
}

// return magnitude of the pseudoscalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr T nrm(PScalar3d<T> ps)
{
    return std::abs(T(ps));
}

// return the reverse of a trivector
template <typename T> inline PScalar3d<T> rev(PScalar3d<T> A)
{
    // the 3d trivector switches sign on reversion
    return PScalar3d<T>(-T(A));
}

// return inverse of the pseudoscalar (A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
// k is the dimension of the space of the pseudoscalar formed by k orthogonal vectors
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> inv(PScalar3d<T> ps)
{
    return -PScalar3d<T>(ps) / sq_nrm(ps);
}

// return the value of the pseudoscalar as value_t (for use in scripting)
template <typename T>
    requires(std::floating_point<T>)
inline constexpr value_t to_val(PScalar3d<T> ps)
{
    return value_t(ps);
}

////////////////////////////////////////////////////////////////////////////////
// PScalar3d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, PScalar3d<T> v)
{
    os << "(" << T(v) << ")";
    return os;
}

} // namespace hd::ga

////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::PScalar3d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::PScalar3d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({})", nested(double(v)));
    }
};