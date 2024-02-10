#ifndef HD_GA_HPP
#define HD_GA_HPP

// author: Daniel Hug, 2024

////////////////////////////////////////////////////////////////////////////////
// provide data types for representing GA in 2D and 3D
// e.g. scalar, vector, bivector, trivector
// and corresponding operations
////////////////////////////////////////////////////////////////////////////////

#include <array>
#include <cstdint> // uint8_t
#include <string>

namespace hd {

using namespace std::literals; // for string initialization

////////////////////////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////////////////////////

// create an algebra with:
// P generators for numbers that square to +1
// N generators for numbers that square to -1
// Z generators for numbers that square to  0
//
// used to provide frequently used values of the algebra
// should be assigned to a const variable
template <uint8_t P, uint8_t N = 0, uint8_t Z = 0>
    requires(P + N + Z >= 2 && P + N + Z <= 4 && P <= 4 && N == 0 &&
             Z == 0) // no implementation for other algebras yet
struct algebra {
    uint8_t p{P};                             // number of +1 generators
    uint8_t n{N};                             // number of -1 generators
    uint8_t z{Z};                             // number of  0 generators
    uint8_t dim_space{P + N + Z};             // dimension of the space
    uint8_t num_components{1 << (P + N + Z)}; // the number of basis components
                                              // == 2^dim_space

    std::array<const uint8_t, P + N + Z + 1> num_components_grade =
        []() -> std::array<const uint8_t, P + N + Z + 1> {
        if constexpr (P + N + Z == 2) return {1, 2, 1};
        if constexpr (P + N + Z == 3) return {1, 3, 3, 1};
        if constexpr (P + N + Z == 4) return {1, 4, 6, 4, 1};
    }();

    std::array<const std::string, 1 << (P + N + Z)> basis_name =
        []() -> std::array<const std::string, (1 << (P + N + Z))> {
        if constexpr (P + N + Z == 2) return {"1"s, "e1"s, "e2"s, "e12"s};
        if constexpr (P + N + Z == 3)
            return {"1"s, "e1"s, "e2"s, "e3"s, "e23"s, "e31"s, "e12"s, "e123"s};
        if constexpr (P + N + Z == 4)
            return {"1"s,   "e1"s,  "e2"s,  "e3"s,   "e0"s,   "e01"s,  "e02"s,  "e03"s,
                    "e12"s, "e31"s, "e23"s, "e032"s, "e013"s, "e021"s, "e123"s, "e0123"s};
    }();
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////


} // namespace hd

#endif // HD_GA_HPP