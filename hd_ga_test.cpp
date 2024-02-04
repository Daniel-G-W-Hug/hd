#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// include functions to be tests
#include "hd_ga.hpp"

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

TEST_SUITE("algebra")
{
    TEST_CASE("2D algebra")
    {
        hd::algebra<2> alg;
        CHECK(alg.p == 2);
        CHECK(alg.n == 0);
        CHECK(alg.z == 0);
        CHECK(alg.dim_space == 2);  // dim == p+n+z
        CHECK(alg.components == 4); // comp == 2^dim
        CHECK(alg.dim_grade.size() == 3);
        fmt::println("2D algebra: dim_grade = {}", fmt::join(alg.dim_grade, ","));
        fmt::println("2D algebra: basis_name = {}", fmt::join(alg.basis_name, ","));
    }

    TEST_CASE("3D algebra")
    {
        hd::algebra<3> alg;
        CHECK(alg.p == 3);
        CHECK(alg.n == 0);
        CHECK(alg.z == 0);
        CHECK(alg.dim_space == 3);  // dim == p+n+z
        CHECK(alg.components == 8); // comp == 2^dim
        CHECK(alg.dim_grade.size() == 4);
        fmt::println("3D algebra: dim_grade = {}", fmt::join(alg.dim_grade, ","));
        fmt::println("3D algebra: basis_name = {}", fmt::join(alg.basis_name, ","));
    }

    TEST_CASE("4D algebra")
    {
        hd::algebra<4> alg;
        CHECK(alg.p == 4);
        CHECK(alg.n == 0);
        CHECK(alg.z == 0);
        CHECK(alg.dim_space == 4);   // dim == p+n+z
        CHECK(alg.components == 16); // comp == 2^dim
        CHECK(alg.dim_grade.size() == 5);
        fmt::println("4D algebra: dim_grade = {}", fmt::join(alg.dim_grade, ","));
        fmt::println("4D algebra: basis_name = {}", fmt::join(alg.basis_name, ","));
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