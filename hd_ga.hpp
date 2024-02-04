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
#include <vector>

namespace hd {

////////////////////////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////////////////////////

const std::vector<std::vector<const uint8_t>> dim_grades{
    {1, 2, 1}, {1, 3, 3, 1}, {1, 4, 6, 4, 1}};

const std::vector<std::vector<std::vector<const std::string>>> basis_names{
    {{"1"}, {"e1", "e2"}, {"e12"}},
    {{"1"}, {"e1", "e2", "e3"}, {"e12", "e23", "e31"}, {"e123"}},
    {{"1"},
     {"e1", "e2", "e3", "e4"},
     {"e12", "e23", "e34", "e41", "e31", "e24"},
     {"e123", "e234", "e341", "e412"},
     {"e1234"}}};

// create an algebra with:
// P generators for numbers that square to +1
// N generators for numbers that square to -1
// Z generators for numbers that square to  0
//
template <uint8_t P, uint8_t N = 0, uint8_t Z = 0>
    requires(P + N + Z >= 2 && P + N + Z <= 4 && P <= 4 && N == 0 &&
             Z == 0) // no implementation for other algebras yet
struct algebra {
    const uint8_t p = P;                         // number of +1 generators
    const uint8_t n = N;                         // number of -1 generators
    const uint8_t z = Z;                         // number of  0 generators
    const uint8_t dim_space = p + n + z;         // dimension of the space
    const uint8_t components = (1 << dim_space); // the number of basis components
                                                 // = 2^dim_space
    const std::vector<const uint8_t> dim_grade{dim_grades[P + N + Z - 2]};
    const std::vector<std::vector<const std::string>> basis_name{
        basis_names[P + N + Z - 2]};
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////


} // namespace hd

#endif // HD_GA_HPP