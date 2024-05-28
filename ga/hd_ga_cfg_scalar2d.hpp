#pragma once

// author: Daniel Hug, 2024

#include <concepts> // std::floating_point

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "hd_ga_cfg_value_t.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Scalar2d<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return the value of the scalar as value_t (for use in scripting)
template <typename T>
    requires(std::floating_point<T>)
inline constexpr value_t to_val(Scalar2d<T> s)
{
    return value_t(s);
}

////////////////////////////////////////////////////////////////////////////////
// Scalar2d<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Scalar2d<T> v)
{
    os << "(" << T(v) << ")";
    return os;
}
} // namespace hd::ga

////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::Scalar2d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Scalar2d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({})", nested(double(v)));
    }
};