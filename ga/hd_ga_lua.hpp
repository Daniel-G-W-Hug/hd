#pragma once

// author: Daniel Hug, 2024

#define SOL_ALL_SAFETIES_ON 1
#include "sol/sol.hpp"

#include <iostream> // cout, cerr

#include "hd_ga.hpp"

////////////////////////////////////////////////////////////////////////////////
// register basic types, geometric operations and constants
// of user-defined types with lua
////////////////////////////////////////////////////////////////////////////////
void register_2d_types(sol::state& lua);
void register_functions(sol::state& lua);
void register_constants(sol::state& lua);


////////////////////////////////////////////////////////////////////////////////
// 2d types
////////////////////////////////////////////////////////////////////////////////
void register_2d_types(sol::state& lua)
{
    using namespace hd::ga;

    lua.new_usertype<scalar2d>(
        "scalar2d",
        sol::constructors<scalar2d(), scalar2d(value_t const&), scalar2d(value_t&&)>());


    lua.new_usertype<pscalar2d>(
        "pscalar2d",
        sol::constructors<pscalar2d(), pscalar2d(value_t const&), pscalar2d(value_t&&)>(),
        sol::meta_function::unary_minus, sol::resolve<pscalar2d(pscalar2d)>(operator-),
        sol::meta_function::addition,
        sol::resolve<pscalar2d(pscalar2d, pscalar2d)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<pscalar2d(pscalar2d, pscalar2d)>(operator-),
        sol::meta_function::multiplication,
        sol::overload( // operator* must be assigned to type with first complex arg
            sol::resolve<mvec2d(pscalar2d, mvec2d const&)>(operator*),
            sol::resolve<mvec2d_e(pscalar2d, mvec2d_e const&)>(operator*),
            sol::resolve<vec2d(pscalar2d, vec2d const&)>(operator*),
            sol::resolve<value_t(pscalar2d, pscalar2d)>(operator*),
            sol::resolve<pscalar2d(pscalar2d, value_t)>(operator*),
            sol::resolve<pscalar2d(value_t, pscalar2d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<pscalar2d(pscalar2d, value_t)>(operator/));


    lua.new_usertype<vec2d>(
        "vec2d", sol::constructors<vec2d(), vec2d(value_t, value_t), vec2d(vec2d)>(), "x",
        &vec2d::x, "y", &vec2d::y, sol::meta_function::unary_minus,
        sol::resolve<vec2d(vec2d const&)>(operator-), sol::meta_function::addition,
        sol::resolve<vec2d(vec2d const&, vec2d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<vec2d(vec2d const&, vec2d const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<vec2d(vec2d const&, value_t)>(operator*),
                      sol::resolve<vec2d(value_t, vec2d const&)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(operator*),
                      sol::resolve<mvec2d(vec2d const&, mvec2d const&)>(operator*),
                      sol::resolve<mvec2d_e(vec2d const&, vec2d const&)>(operator*),
                      sol::resolve<vec2d(vec2d const&, pscalar2d)>(operator*)),
        sol::meta_function::division,
        sol::resolve<vec2d(vec2d const&, value_t)>(operator/));


    lua.new_usertype<mvec2d_e>(
        "mvec2d_e",
        sol::constructors<mvec2d_e(), mvec2d_e(value_t, value_t), mvec2d_e(scalar2d),
                          mvec2d_e(pscalar2d), mvec2d_e(scalar2d, pscalar2d),
                          mvec2d_e(mvec2d_e)>(),
        "c0", &mvec2d_e::c0, "c1", &mvec2d_e::c1, sol::meta_function::unary_minus,
        sol::resolve<mvec2d_e(mvec2d_e const&)>(operator-), sol::meta_function::addition,
        sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&, value_t)>(operator*),
                      sol::resolve<mvec2d_e(value_t, mvec2d_e const&)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<vec2d(mvec2d_e const&, vec2d const&)>(operator*),
                      sol::resolve<mvec2d_e(mvec2d_e const&, mvec2d_e const&)>(operator*),
                      sol::resolve<mvec2d_e(mvec2d_e const&, pscalar2d)>(operator*),
                      sol::resolve<mvec2d(mvec2d_e const&, mvec2d const&)>(operator*),
                      sol::resolve<vec2d(mvec2d_e const&, vec2d const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec2d_e(mvec2d_e const&, value_t)>(operator/));


    lua.new_usertype<mvec2d>(
        "mvec2d",
        sol::constructors<mvec2d(), mvec2d(value_t, value_t, value_t, value_t),
                          mvec2d(scalar2d), mvec2d(vec2d), mvec2d(pscalar2d),
                          mvec2d(scalar2d, pscalar2d), mvec2d(mvec2d_e),
                          mvec2d(mvec2d)>(),
        "c0", &mvec2d::c0, "c1", &mvec2d::c1, "c2", &mvec2d::c2, "c3", &mvec2d::c3,
        sol::meta_function::unary_minus, sol::resolve<mvec2d(mvec2d const&)>(operator-),
        sol::meta_function::addition,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator+),
        sol::meta_function::subtraction,
        sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator-),
        sol::meta_function::multiplication,
        sol::overload(sol::resolve<mvec2d(mvec2d const&, value_t)>(operator*),
                      sol::resolve<mvec2d(value_t, mvec2d const&)>(operator*),
                      // operator* must be assigned to type with first complex arg
                      sol::resolve<mvec2d(mvec2d const&, vec2d const&)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d const&)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, pscalar2d)>(operator*),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d_e const&)>(operator*)),
        sol::meta_function::division,
        sol::resolve<mvec2d(mvec2d const&, value_t)>(operator/));
}

////////////////////////////////////////////////////////////////////////////////
// functions commonly used across 2d and 3d types
////////////////////////////////////////////////////////////////////////////////
void register_functions(sol::state& lua)
{
    using namespace hd::ga;

    lua.set_function("dot", sol::resolve<value_t(vec2d const&, vec2d const&)>(dot));

    lua.set_function("sq_nrm",
                     sol::overload(sol::resolve<value_t(vec2d const&)>(sq_nrm),
                                   sol::resolve<value_t(mvec2d_e const&)>(sq_nrm),
                                   sol::resolve<value_t(mvec2d const&)>(sq_nrm),
                                   sol::resolve<value_t(pscalar2d)>(sq_nrm)));

    lua.set_function("nrm", sol::overload(sol::resolve<value_t(vec2d const&)>(nrm),
                                          sol::resolve<value_t(mvec2d_e const&)>(nrm),
                                          sol::resolve<value_t(mvec2d const&)>(nrm),
                                          sol::resolve<value_t(pscalar2d)>(nrm)));

    lua.set_function("rev", sol::overload(sol::resolve<mvec2d_e(mvec2d_e const&)>(rev),
                                          sol::resolve<mvec2d(mvec2d const&)>(rev)));

    lua.set_function("conj", sol::overload(sol::resolve<mvec2d(mvec2d const&)>(conj)));

    lua.set_function("unitized",
                     sol::overload(sol::resolve<vec2d(vec2d const&)>(unitized),
                                   sol::resolve<mvec2d_e(mvec2d_e const&)>(unitized),
                                   sol::resolve<mvec2d(mvec2d const&)>(unitized)));

    lua.set_function("inv", sol::overload(sol::resolve<vec2d(vec2d const&)>(inv),
                                          sol::resolve<mvec2d_e(mvec2d_e const&)>(inv),
                                          sol::resolve<mvec2d(mvec2d const&)>(inv),
                                          sol::resolve<pscalar2d(pscalar2d)>(inv)));

    lua.set_function("wdg", sol::resolve<pscalar2d(vec2d const&, vec2d const&)>(wdg));

    ////////////////////////////////////////////////////////////////////////////////
    // projections, rejections and reflections
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("project_onto",
                     sol::resolve<vec2d(vec2d const&, vec2d const&)>(project_onto));

    lua.set_function(
        "project_onto_unitized",
        sol::resolve<vec2d(vec2d const&, vec2d const&)>(project_onto_unitized));

    lua.set_function("reject_from",
                     sol::resolve<vec2d(vec2d const&, vec2d const&)>(reject_from));

    lua.set_function(
        "reject_from_unitized",
        sol::resolve<vec2d(vec2d const&, vec2d const&)>(reject_from_unitized));

    lua.set_function("reflect_on_hyp",
                     sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on_hyp));

    lua.set_function("reflect_on_vec",
                     sol::resolve<vec2d(vec2d const&, vec2d const&)>(reflect_on_vec));

    ////////////////////////////////////////////////////////////////////////////////
    // angles and rotations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("angle", sol::resolve<value_t(vec2d const&, vec2d const&)>(angle));

    lua.set_function("angle_to_re", sol::resolve<value_t(mvec2d_e const&)>(angle_to_re));

    lua.set_function("exp", sol::resolve<mvec2d_e(pscalar2d, value_t)>(exp));

    lua.set_function("rotor", sol::resolve<mvec2d_e(pscalar2d, value_t)>(rotor));

    lua.set_function(
        "rotate",
        sol::overload(sol::resolve<vec2d(vec2d const&, mvec2d_e const&)>(rotate),
                      sol::resolve<mvec2d(mvec2d const&, mvec2d_e const&)>(rotate)));

    ////////////////////////////////////////////////////////////////////////////////
    // dualization operations
    ////////////////////////////////////////////////////////////////////////////////

    lua.set_function("dual", sol::overload(sol::resolve<scalar2d(pscalar2d)>(dual),
                                           sol::resolve<pscalar2d(scalar2d)>(dual),
                                           sol::resolve<vec2d(vec2d const&)>(dual),
                                           sol::resolve<mvec2d_e(mvec2d_e const&)>(dual),
                                           sol::resolve<mvec2d(mvec2d const&)>(dual)));

    ////////////////////////////////////////////////////////////////////////////////
    // grade operations on multivectors
    ////////////////////////////////////////////////////////////////////////////////
    lua.set_function("gr0", sol::overload(sol::resolve<scalar2d(mvec2d_e const&)>(gr0),
                                          sol::resolve<scalar2d(mvec2d const&)>(gr0)));

    lua.set_function("gr1", sol::resolve<vec2d(mvec2d const&)>(gr1));

    lua.set_function("gr2", sol::overload(sol::resolve<pscalar2d(mvec2d_e const&)>(gr2),
                                          sol::resolve<pscalar2d(mvec2d const&)>(gr2)));

    ////////////////////////////////////////////////////////////////////////////////
    // common helper functions for scripting in lua
    ////////////////////////////////////////////////////////////////////////////////

    // convert scalars & pscalars into numeric values for further calculations
    lua.set_function("to_val", sol::overload(sol::resolve<value_t(scalar2d)>(to_val),
                                             sol::resolve<value_t(pscalar2d)>(to_val)));

    lua.set_function("rad_to_deg", &rad_to_deg);
    lua.set_function("deg_to_rad", &deg_to_rad);
}

////////////////////////////////////////////////////////////////////////////////
// make defined constants available as global variables in lua
////////////////////////////////////////////////////////////////////////////////
void register_constants(sol::state& lua)
{
    using namespace hd::ga;

    // 2D constants
    lua["eps"] = eps;
    lua["e1_2d"] = e1_2d;
    lua["e2_2d"] = e2_2d;
    lua["e1m_2d"] = e1m_2d;
    lua["e2m_2d"] = e2m_2d;
    lua["I_2d"] = I_2d;
    lua["Im_2d"] = Im_2d;
    lua["Im_2d_E"] = Im_2d_E;
}