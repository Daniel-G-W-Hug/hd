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

// definition of scalar and pscalar located here, because other
// modules like vec2d, vec3d, 2d_ops, 3d_ops etc. just need the type defintion,
// but not the functionality itself (reduce dependencies between files)

struct scalar2d_tag {};
struct scalar3d_tag {};
struct pscalar2d_tag {};
struct pscalar3d_tag {};

template <typename T> using Scalar2d = Strong_t<T, scalar2d_tag>;
template <typename T> using Scalar3d = Strong_t<T, scalar3d_tag>;
template <typename T> using PScalar2d = Strong_t<T, pscalar2d_tag>;
template <typename T> using PScalar3d = Strong_t<T, pscalar3d_tag>;

} // namespace hd::ga