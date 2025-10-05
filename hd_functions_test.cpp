#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// include functions to be tests
#include "hd_functions.hpp"

using namespace hd; // find all functions in hd:: namespace

TEST_SUITE("fact(n):")
{
    TEST_CASE("fact(n): specific values")
    {
        CHECK(fact(0) == 1.0);
        CHECK(fact(1) == 1.0);
        CHECK(fact(2) == 2.0);
        CHECK(fact(3) == 6.0);
        CHECK(fact(10) == 3628800.);
    }

    TEST_CASE("fact(n): throws on n<0") { CHECK_THROWS(fact(-1)); }

    TEST_CASE(("Kronecker-Delta and Levi-Cevita permutation symbol eps(...)"))
    {
        CHECK(kronecker(0, 0) == 1);
        CHECK(kronecker(1, 0) == 0);
        CHECK(kronecker(0, 1) == 0);
        CHECK(kronecker(1, 1) == 1);

        CHECK(kronecker<>(0, 0) == 1);
        CHECK(kronecker<>(1, 0) == 0);
        CHECK(kronecker<>(0, 1) == 0);
        CHECK(kronecker<>(1, 1) == 1);

        CHECK(kronecker<int>(0, 0) == 1);
        CHECK(kronecker<int>(1, 0) == 0);
        CHECK(kronecker<int>(0, 1) == 0);
        CHECK(kronecker<int>(1, 1) == 1);

        CHECK(kronecker<double>(0, 0) == 1.0);
        CHECK(kronecker<double>(1, 0) == 0.0);
        CHECK(kronecker<double>(0, 1) == 0.0);
        CHECK(kronecker<double>(1, 1) == 1.0);

        CHECK(kronecker<float>(0, 0) == doctest::Approx(1.0f));
        CHECK(kronecker<float>(1, 0) == doctest::Approx(0.0f));
        CHECK(kronecker<float>(0, 1) == doctest::Approx(0.0f));
        CHECK(kronecker<float>(1, 1) == doctest::Approx(1.0f));

        // 3D cases:

        // even permutations
        CHECK(eps(0, 1, 2) == 1);
        CHECK(eps(1, 2, 0) == 1);
        CHECK(eps(2, 0, 1) == 1);

        // odd permutations
        CHECK(eps(1, 0, 2) == -1);
        CHECK(eps(2, 1, 0) == -1);
        CHECK(eps(0, 2, 1) == -1);

        // double indices
        CHECK(eps(0, 0, 1) == 0);
        CHECK(eps(0, 1, 1) == 0);
        CHECK(eps(1, 1, 2) == 0);

        // 4D cases:
        CHECK(eps(0, 1, 2, 3) == 1);

        // simple transposition (odd
        CHECK(eps(1, 0, 2, 3) == -1);
        CHECK(eps(0, 2, 1, 3) == -1);
        CHECK(eps(0, 1, 3, 2) == -1);

        // double transposition (even)
        CHECK(eps(1, 0, 3, 2) == 1);


        // double indices
        CHECK(eps(0, 0, 1, 2) == 0);
        CHECK(eps(0, 1, 2, 2) == 0);

        // 5D case:
        CHECK(eps(0, 1, 2, 3, 4) == 1);
    }
}

TEST_CASE("Levi-Civita Symbol in 2D")
{
    SUBCASE("Even permutation (identity)") { CHECK(eps(0, 1) == 1); }

    SUBCASE("Odd permutation (one transposition)") { CHECK(eps(1, 0) == -1); }

    SUBCASE("Repeated indices")
    {
        CHECK(eps(0, 0) == 0);
        CHECK(eps(1, 1) == 0);
    }
}

TEST_CASE("Levi-Civita Symbol in 3D - All permutations")
{
    SUBCASE("Even permutations (+1)")
    {
        CHECK(eps(0, 1, 2) == 1); // Identity
        CHECK(eps(1, 2, 0) == 1); // Cyclic rotation
        CHECK(eps(2, 0, 1) == 1); // Cyclic rotation
    }

    SUBCASE("Odd permutations (-1)")
    {
        CHECK(eps(0, 2, 1) == -1); // One transposition
        CHECK(eps(1, 0, 2) == -1); // One transposition
        CHECK(eps(2, 1, 0) == -1); // One transposition
    }

    SUBCASE("Repeated indices (0)")
    {
        CHECK(eps(0, 0, 1) == 0);
        CHECK(eps(0, 1, 1) == 0);
        CHECK(eps(1, 1, 2) == 0);
        CHECK(eps(0, 0, 0) == 0);
        CHECK(eps(1, 0, 1) == 0);
        CHECK(eps(2, 2, 0) == 0);
    }
}

TEST_CASE("Levi-Civita Symbol in 3D - Alternative index sets")
{
    SUBCASE("Using indices {0,1,2}")
    {
        CHECK(eps(0, 1, 2) == 1);
        CHECK(eps(2, 0, 1) == 1);
        CHECK(eps(1, 0, 2) == -1);
    }

    SUBCASE("Using indices {1,2,3}")
    {
        CHECK(eps(1, 2, 3) == 1);
        CHECK(eps(2, 3, 1) == 1);
        CHECK(eps(3, 1, 2) == 1);
        CHECK(eps(1, 3, 2) == -1);
    }

    SUBCASE("Invalid indices should throw")
    {
        CHECK_THROWS(eps(0, 2, 4)); // Nicht aufeinanderfolgend
        CHECK_THROWS(eps(0, 1, 5)); // Nicht aufeinanderfolgend
        CHECK_THROWS(eps(2, 3, 5)); // Nicht aufeinanderfolgend
    }
}

TEST_CASE("Levi-Civita Symbol - Antisymmetry property")
{
    SUBCASE("2D: Swapping two indices changes sign") { CHECK(eps(0, 1) == -eps(1, 0)); }

    SUBCASE("3D: Swapping two indices changes sign")
    {
        CHECK(eps(0, 1, 2) == -eps(1, 0, 2));
        CHECK(eps(0, 1, 2) == -eps(0, 2, 1));
        CHECK(eps(0, 1, 2) == -eps(2, 1, 0));
        CHECK(eps(1, 2, 0) == -eps(2, 1, 0));
    }
}

TEST_CASE("Levi-Civita Symbol - Type compatibility")
{
    SUBCASE("2D with size_t")
    {
        CHECK(eps(size_t(0), size_t(1)) == 1);
        CHECK(eps(size_t(1), size_t(0)) == -1);
    }

    SUBCASE("3D with size_t")
    {
        CHECK(eps(size_t(0), size_t(1), size_t(2)) == 1);
        CHECK(eps(size_t(2), size_t(1), size_t(0)) == -1);
    }

    SUBCASE("3D with mixed int and size_t")
    {
        CHECK(eps(0, size_t(1), 2) == 1);
        CHECK(eps(size_t(1), 0, size_t(2)) == -1);
    }
}
