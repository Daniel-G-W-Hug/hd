#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// include functions to be tested
#include "hd_determinant.hpp"

using namespace hd; // find all functions in hd:: namespace

TEST_SUITE("Determinant calculations:")
{
    TEST_CASE("det() with vector<vector<T>>: 2x2 matrices")
    {
        SUBCASE("Identity matrix")
        {
            std::vector<std::vector<double>> I = {{1.0, 0.0}, {0.0, 1.0}};
            CHECK(det(I) == doctest::Approx(1.0));
        }

        SUBCASE("Simple 2x2 matrix")
        {
            std::vector<std::vector<double>> A = {{1.0, 2.0}, {3.0, 4.0}};
            // det = 1*4 - 2*3 = -2
            CHECK(det(A) == doctest::Approx(-2.0));
        }

        SUBCASE("Integer 2x2 matrix")
        {
            std::vector<std::vector<int>> A = {{2, 3}, {1, 4}};
            // det = 2*4 - 3*1 = 5
            CHECK(det(A) == 5);
        }

        SUBCASE("Singular 2x2 matrix")
        {
            std::vector<std::vector<double>> A = {{1.0, 2.0}, {2.0, 4.0}};
            CHECK(det(A) == doctest::Approx(0.0));
        }
    }

    TEST_CASE("det() with vector<vector<T>>: 3x3 matrices")
    {
        SUBCASE("Identity matrix")
        {
            std::vector<std::vector<double>> I = {
                {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
            CHECK(det(I) == doctest::Approx(1.0));
        }

        SUBCASE("Simple 3x3 matrix")
        {
            std::vector<std::vector<double>> A = {
                {1.0, 2.0, 3.0}, {0.0, 1.0, 4.0}, {5.0, 6.0, 0.0}};
            // det = 1*(1*0 - 4*6) - 2*(0*0 - 4*5) + 3*(0*6 - 1*5)
            //     = 1*(-24) - 2*(-20) + 3*(-5)
            //     = -24 + 40 - 15 = 1
            CHECK(det(A) == doctest::Approx(1.0));
        }

        SUBCASE("Integer 3x3 matrix from user example")
        {
            std::vector<std::vector<int>> A = {{0, 2, 6}, {1, 8, 4}, {5, 2, 7}};
            // Using cofactor expansion along first row:
            // det = 0*(8*7-4*2) - 2*(1*7-4*5) + 6*(1*2-8*5)
            //     = 0 - 2*(7-20) + 6*(2-40)
            //     = 0 - 2*(-13) + 6*(-38)
            //     = 0 + 26 - 228 = -202
            CHECK(det(A) == -202);
        }

        SUBCASE("Singular 3x3 matrix")
        {
            std::vector<std::vector<double>> A = {
                {1.0, 2.0, 3.0}, {2.0, 4.0, 6.0}, {1.0, 1.0, 1.0}};
            CHECK(det(A) == doctest::Approx(0.0));
        }
    }

    TEST_CASE("det() with vector<vector<T>>: 4x4 matrices")
    {
        SUBCASE("Identity matrix")
        {
            std::vector<std::vector<double>> I = {{1.0, 0.0, 0.0, 0.0},
                                                  {0.0, 1.0, 0.0, 0.0},
                                                  {0.0, 0.0, 1.0, 0.0},
                                                  {0.0, 0.0, 0.0, 1.0}};
            CHECK(det(I) == doctest::Approx(1.0));
        }

        SUBCASE("Diagonal matrix")
        {
            std::vector<std::vector<double>> D = {{2.0, 0.0, 0.0, 0.0},
                                                  {0.0, 3.0, 0.0, 0.0},
                                                  {0.0, 0.0, 4.0, 0.0},
                                                  {0.0, 0.0, 0.0, 5.0}};
            // det = 2*3*4*5 = 120
            CHECK(det(D) == doctest::Approx(120.0));
        }

        SUBCASE("General 4x4 matrix")
        {
            std::vector<std::vector<double>> A = {{1.0, 2.0, 0.0, 1.0},
                                                  {3.0, 1.0, 2.0, 0.0},
                                                  {0.0, 1.0, 1.0, 2.0},
                                                  {2.0, 0.0, 1.0, 1.0}};
            // det = -14
            CHECK(det(A) == doctest::Approx(-14.0));
        }

        SUBCASE("Singular 4x4 matrix")
        {
            std::vector<std::vector<double>> A = {{1.0, 2.0, 3.0, 4.0},
                                                  {2.0, 4.0, 6.0, 8.0},
                                                  {1.0, 1.0, 1.0, 1.0},
                                                  {0.0, 1.0, 2.0, 3.0}};
            CHECK(det(A) == doctest::Approx(0.0));
        }
    }

    TEST_CASE("det() with mdspan: 2x2 matrices")
    {
        SUBCASE("Identity matrix")
        {
            double data[] = {1.0, 0.0, 0.0, 1.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> I(data, 2, 2);
            CHECK(det<double>(I) == doctest::Approx(1.0));
        }

        SUBCASE("Simple 2x2 matrix")
        {
            double data[] = {1.0, 2.0, 3.0, 4.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> A(data, 2, 2);
            CHECK(det<double>(A) == doctest::Approx(-2.0));
        }

        SUBCASE("Singular 2x2 matrix")
        {
            double data[] = {1.0, 2.0, 2.0, 4.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> A(data, 2, 2);
            CHECK(det<double>(A) == doctest::Approx(0.0));
        }
    }

    TEST_CASE("det() with mdspan: 3x3 matrices")
    {
        SUBCASE("Identity matrix")
        {
            double data[] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> I(data, 3, 3);
            CHECK(det<double>(I) == doctest::Approx(1.0));
        }

        SUBCASE("Simple 3x3 matrix")
        {
            double data[] = {1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> A(data, 3, 3);
            CHECK(det<double>(A) == doctest::Approx(1.0));
        }

        SUBCASE("Integer 3x3 matrix from user example")
        {
            double data[] = {0.0, 2.0, 6.0, 1.0, 8.0, 4.0, 5.0, 2.0, 7.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> A(data, 3, 3);
            CHECK(det<double>(A) == doctest::Approx(-202.0));
        }
    }

    TEST_CASE("det() with mdspan: 4x4 matrices")
    {
        SUBCASE("Identity matrix")
        {
            double data[] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                             0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> I(data, 4, 4);
            CHECK(det<double>(I) == doctest::Approx(1.0));
        }

        SUBCASE("Diagonal matrix")
        {
            double data[] = {2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0,
                             0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 5.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> D(data, 4, 4);
            CHECK(det<double>(D) == doctest::Approx(120.0));
        }

        SUBCASE("General 4x4 matrix")
        {
            double data[] = {1.0, 2.0, 0.0, 1.0, 3.0, 1.0, 2.0, 0.0,
                             0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 1.0, 1.0};
            Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> A(data, 4, 4);
            CHECK(det<double>(A) == doctest::Approx(-14.0));
        }
    }

    TEST_CASE("det() with vector<vector<T>>: 5x5 matrix")
    {
        SUBCASE("Random 5x5 matrix with entries in [-1, 1]")
        {
            std::vector<std::vector<double>> A = {
                {-0.250920, 0.901429, 0.463988, 0.197317, -0.687963},
                {-0.688011, -0.883833, 0.732352, 0.202230, 0.416145},
                {-0.958831, 0.939820, 0.664885, -0.575322, -0.636350},
                {-0.633191, -0.391516, 0.049513, -0.136110, -0.417542},
                {0.223706, -0.721012, -0.415711, -0.267276, -0.087860}};
            // Verified with numpy: det = 0.3025895007082695
            CHECK(det(A) == doctest::Approx(0.302589501).epsilon(1e-6));
        }
    }

    TEST_CASE("det() error handling")
    {
        SUBCASE("Empty matrix throws")
        {
            std::vector<std::vector<double>> empty;
            CHECK_THROWS(det(empty));
        }

        SUBCASE("Non-square matrix throws")
        {
            std::vector<std::vector<double>> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
            CHECK_THROWS(det(A));
        }
    }
}
