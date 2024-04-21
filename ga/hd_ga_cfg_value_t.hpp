#pragma once

// author: Daniel Hug, 2024

#include "hd_ga_cfg_strong_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions for easy use
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

// select the floating point type used for scalars, vector and bivector components
using value_t = float;
// using value_t = double;

struct scalar_tag {};
struct pscalar2d_tag {};
struct pscalar3d_tag {};

template <typename T> using Scalar = Strong_t<T, scalar_tag>;
template <typename T> using PScalar2d = Strong_t<T, pscalar2d_tag>;
template <typename T> using PScalar3d = Strong_t<T, pscalar3d_tag>;

} // namespace hd::ga