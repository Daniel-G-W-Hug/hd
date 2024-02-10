#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// include functions to be tests
#include "hd_ga.hpp"

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

TEST_SUITE("algebra")
{
    TEST_CASE("2d_ega")
    {
        // 2d euklidean geometric algebra
        const hd::algebra<2> alg;
        CHECK(alg.p == 2);
        CHECK(alg.n == 0);
        CHECK(alg.z == 0);
        CHECK(alg.dim_space == 2);                   // dim_space == p+n+z
        CHECK(alg.num_components == 4);              // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 3); // == dim_space + 1
        fmt::println("2d_ega: dim_grade = {}", fmt::join(alg.num_components_grade, ", "));
        fmt::println("2d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("3d_ega")
    {
        // 3d euklidean geometric algebra
        const hd::algebra<3> alg;
        CHECK(alg.p == 3);
        CHECK(alg.n == 0);
        CHECK(alg.z == 0);
        CHECK(alg.dim_space == 3);                   // dim_space == p+n+z
        CHECK(alg.num_components == 8);              // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("3d_ega: dim_grade = {}", fmt::join(alg.num_components_grade, ", "));
        fmt::println("3d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("4d_ega")
    {
        // 4d euklidean geometric algebra
        const hd::algebra<4> alg;
        CHECK(alg.p == 4);
        CHECK(alg.n == 0);
        CHECK(alg.z == 0);
        CHECK(alg.dim_space == 4);                   // dim_space == p+n+z
        CHECK(alg.num_components == 16);             // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("4d_ega: dim_grade = {}", fmt::join(alg.num_components_grade, ", "));
        fmt::println("4d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }
}


// TEST_SUITE("fact(n):")
// {
//     TEST_CASE("fact(n): specific values")
//     {
//         CHECK(hd::fact(0) == 1.0);
//         CHECK(hd::fact(1) == 1.0);
//         CHECK(hd::fact(2) == 2.0);
//         CHECK(hd::fact(3) == 6.0);
//         CHECK(hd::fact(10) == 3628800.);
//     }
//     TEST_CASE("fact(n): throws on n<0")
//     {
//         CHECK_THROWS(hd::fact(-1));
//     }
// }