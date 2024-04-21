#pragma once

#include "hd_ga_cfg_bivec3d.hpp" // BiVec3d<T>
#include "hd_ga_cfg_mvec2d.hpp"  // Vec2d<T>
#include "hd_ga_cfg_mvec3d.hpp"  // Vec3d<T>
#include "hd_ga_cfg_strong_t.hpp"
#include "hd_ga_cfg_value_t.hpp" // defines value_t either as float oder double
#include "hd_ga_cfg_vec2d.hpp"   // Vec2d<T>
#include "hd_ga_cfg_vec3d.hpp"   // Vec3d<T>

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// 2D constants
////////////////////////////////////////////////////////////////////////////////

auto e1_2d = Vec2d<value_t>{1.0, 0.0};
auto e2_2d = Vec2d<value_t>{0.0, 1.0};
auto e1m_2d = MVec2d<value_t>{e1_2d};
auto e2m_2d = MVec2d<value_t>{e2_2d};

auto I_2d = PScalar2d<value_t>(1.0);
auto Im_2d = MVec2d<value_t>{I_2d};

////////////////////////////////////////////////////////////////////////////////
// 3D constants
////////////////////////////////////////////////////////////////////////////////

auto e1_3d = Vec3d<value_t>{1.0, 0.0, 0.0};
auto e2_3d = Vec3d<value_t>{0.0, 1.0, 0.0};
auto e3_3d = Vec3d<value_t>{0.0, 0.0, 1.0};
auto e1m_3d = MVec3d<value_t>{e1_3d};
auto e2m_3d = MVec3d<value_t>{e2_3d};
auto e3m_3d = MVec3d<value_t>{e3_3d};

auto e23_3d = BiVec3d<value_t>{1.0, 0.0, 0.0};
auto e31_3d = BiVec3d<value_t>{0.0, 1.0, 0.0};
auto e12_3d = BiVec3d<value_t>{0.0, 0.0, 1.0};
auto e23m_3d = MVec3d<value_t>{e23_3d};
auto e31m_3d = MVec3d<value_t>{e31_3d};
auto e12m_3d = MVec3d<value_t>{e12_3d};

auto I_3d = PScalar3d<value_t>(1.0);
auto Im_3d = MVec3d<value_t>{I_3d};

} // namespace hd::ga