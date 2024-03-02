#pragma once

// author: Daniel Hug, 2024

// this is a convenience header to include all required headers consistently

////////////////////////////////////////////////////////////////////////////////
// provide data types for representing GA in 2D and 3D
// e.g. scalar, vector, bivector, etc.
// and corresponding operations
////////////////////////////////////////////////////////////////////////////////

#include "hd_ga_cfg_algebra.hpp" // algebra
#include "hd_ga_cfg_value_t.hpp" // default type for scalars, vector and bivector components
#include "hd_ga_cfg_vec2d.hpp"   // vec2d<T>
#include "hd_ga_usr_types.hpp"   // consistent user types for all compontents


// struct Frame {
//     Frame(const Vector2d& v1, const Vector2d& v2) : e1(v1), e2(v2) {}
//     Vector2d e1;
//     Vector2d e2;
// };

// T& s{&c[0]};   // scalar part
// T& x{&c[1]};   // vector part, basis vector e1
// T& y{&c[2]};   // vector part, basis vector e2
// T& ps{&c[3]}; // pseudo scalar part