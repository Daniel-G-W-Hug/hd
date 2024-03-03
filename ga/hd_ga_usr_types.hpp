#pragma once

#include "hd_ga_cfg_value_t.hpp"
#include "hd_ga_cfg_vec2d.hpp" // vec2d<T>

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions for easy use
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

// scalar types
using scalar_t = value_t;
using pseudoscalar_t = value_t;

// vector and bivector types
using vec2d = Vec2d<value_t>;

value_t eps{2.0 * std::numeric_limits<value_t>::epsilon()};

} // namespace hd::ga