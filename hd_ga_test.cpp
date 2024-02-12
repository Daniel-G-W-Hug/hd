#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// include functions to be tests
#include "hd_ga.hpp"

#include "fmt/format.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

TEST_SUITE("algebra")
{

    using namespace hd;

    TEST_CASE("2d_ega")
    {
        // 2d euklidean geometric algebra
        const Algebra<2> alg;
        CHECK(alg.p() == 2);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 2);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 4);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 3); // == dim_space + 1
        fmt::println("2d_ega: dim_grade = {}", fmt::join(alg.num_components_grade, ", "));
        fmt::println("2d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("3d_ega")
    {
        // 3d euklidean geometric algebra
        const Algebra<3> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 3);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 8);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("3d_ega: dim_grade = {}", fmt::join(alg.num_components_grade, ", "));
        fmt::println("3d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("4d_ega")
    {
        // 4d euklidean geometric algebra
        const Algebra<4> alg;
        CHECK(alg.p() == 4);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("4d_ega: dim_grade = {}", fmt::join(alg.num_components_grade, ", "));
        fmt::println("4d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("Vec2d default init")
    {
        // default initialization
        Vec2d v;
        CHECK(v.x == 0.0f);
        CHECK(v.y == 0.0f);
    }
    TEST_CASE("Vec2d with curly braced intializer")
    {
        // default initialization
        Vec2d v{0.0f, 0.0f};
        CHECK(v.x == 0.0f);
        CHECK(v.y == 0.0f);
    }
    TEST_CASE("Vec2d cp ctor & cp assign incl. type deduction")
    {
        // default initialization
        Vec2d v1{1.0, 2.0}; // init with double (type deduction)
        Vec2d v2{v1};       // cp ctor
        Vec2d v3 = v2;      // cp assign
        CHECK(v1.x == 1.0);
        CHECK(v1.y == 2.0);
        CHECK(v2.x == 1.0);
        CHECK(v2.y == 2.0);
        CHECK(v3.x == 1.0);
        CHECK(v3.y == 2.0);
    }

    TEST_CASE("Vec2d comparison")
    {
        // default initialization
        Vector2d v1{1.0, 2.0};
        Vector2d v2{2.0, 4.0};
        Vector2d v3{1.0, -2.0000001};
        Vector2d v4{v1};

        // fmt::println("v1 = {}", fmt::streamed(v1));
        // fmt::println("v2 = {}", fmt::streamed(v2));
        // fmt::println("v3 = {}", fmt::streamed(v3));
        // fmt::println("v4 = {}", fmt::streamed(v4));

        // fmt::println("SquaredNorm(v1) = {:e}", SquaredNorm(v1));
        // fmt::println("SquaredNorm(v2) = {:e}", SquaredNorm(v2));
        // fmt::println("SquaredNorm(v3) = {:e}", SquaredNorm(v3));
        // fmt::println("SquaredNorm(v4) = {:e}", SquaredNorm(v4));
        // fmt::println("SquaredNorm(Delta) = {:.20f}", SquaredNorm(v1) -
        // SquaredNorm(v3));

        CHECK(v1 == v4); // comparison (equality)
        CHECK(v1 != v2); // comparison (inequality)
        CHECK(v1 < v2);  // comparison (less than)
        CHECK(v2 >= v1); // comparison (greater than or equal)
        CHECK(v3 == v1); // comparison (eqality)
    }

    TEST_CASE("Vec2d operations")
    {
        Vector2d v1{2.0, 1.0};
        Vector2d v2 = Normalize(v1);

        fmt::println("v1 = {}, Norm(v1) = {}", fmt::streamed(v1), Norm(v1));
        fmt::println("v2 = Normalize(v1) = {}, Norm(v2) = {}", fmt::streamed(v2),
                     Norm(v2));

        CHECK(std::abs(SquaredNorm(v1) - 5.0) < 1.e-7); // comparison (equality)
        CHECK(std::abs(SquaredNorm(v2) - 1.0) < 1.e-7); // comparison (inequality)
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