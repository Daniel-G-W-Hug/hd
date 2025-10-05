#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "fmt/format.h"
#include "fmt/ranges.h"

#include "hd_solver.hpp"

// Test cases generated with generate_test_solver.py using NumPy as reference

double const eps{1.e-15};

TEST_SUITE("LU Solver Tests:")
{
    TEST_CASE("LU solver: original simple system")
    {
        // storage for matrix, rhs & permutation
        std::array m_s{1., 2., 3., 0., 4., 1., 0., 0., 1.};
        std::array rhs_s{1., 1., 1.};
        std::array<int, 3> m_perm_s;

        // setup the corresponding mdarray views onto the data
        auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
        auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
        auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

        // LU decomposition of matrix
        hd::lu_decomp(m, m_perm);

        // solution by backsubstition of rhs => solution is returned in rhs
        hd::lu_backsubs(m, m_perm, rhs);

        for (size_t i = 0; i < m.extent(0); ++i) {
            for (size_t j = 0; j < m.extent(1); ++j) {
                fmt::print("m[{},{}] == {}, ", i, j, m[i, j]);
            }
            fmt::print("\n");
        }

        for (size_t i = 0; i < m.extent(0); ++i)
            fmt::print("rhs[{}] == {},", i, rhs[i]);
        fmt::print("\n");

        CHECK(std::abs(rhs[0] - (-2.0)) < eps);
        CHECK(std::abs(rhs[1] - 0.0) < eps);
        CHECK(std::abs(rhs[2] - 1.0) < eps);
    }

    TEST_CASE("LU solver: basic systems")
    {
        SUBCASE("Simple 2x2 system")
        {
            // Matrix A
            std::array m_s{2.0000000000000000e+00, 1.0000000000000000e+00,
                           1.0000000000000000e+00, 3.0000000000000000e+00};
            // Right-hand side b
            std::array rhs_s{5.0000000000000000e+00, 6.0000000000000000e+00};
            // Expected solution x
            std::array x_expected{1.8000000000000000e+00, 1.3999999999999999e+00};
            std::array<int, 2> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 2, 2>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 2>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 2>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 2; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }

        SUBCASE("Upper triangular 3x3 system")
        {
            // Matrix A
            std::array m_s{
                1.0000000000000000e+00, 2.0000000000000000e+00, 3.0000000000000000e+00,
                0.0000000000000000e+00, 4.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00};
            // Right-hand side b
            std::array rhs_s{1.0000000000000000e+00, 1.0000000000000000e+00,
                             1.0000000000000000e+00};
            // Expected solution x
            std::array x_expected{-2.0000000000000000e+00, 0.0000000000000000e+00,
                                  1.0000000000000000e+00};
            std::array<int, 3> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 3; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }

        SUBCASE("General 4x4 system")
        {
            // Condition number: 1.042765e+02
            // Matrix A
            std::array m_s{
                2.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 4.0000000000000000e+00, 3.0000000000000000e+00,
                3.0000000000000000e+00, 1.0000000000000000e+00, 8.0000000000000000e+00,
                7.0000000000000000e+00, 9.0000000000000000e+00, 5.0000000000000000e+00,
                6.0000000000000000e+00, 7.0000000000000000e+00, 9.0000000000000000e+00,
                8.0000000000000000e+00};
            // Right-hand side b
            std::array rhs_s{1.0000000000000000e+00, 2.0000000000000000e+00,
                             3.0000000000000000e+00, 4.0000000000000000e+00};
            // Expected solution x: [1.0, 0.5, -1.5, 1.0]
            std::array x_expected{1.0000000000000007e+00, 4.9999999999999922e-01,
                                  -1.5000000000000002e+00, 1.0000000000000004e+00};
            std::array<int, 4> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 4, 4>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 4>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 4>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 4; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }

        SUBCASE("Random 5x5 diagonally dominant system")
        {
            // Matrix A
            std::array m_s{
                1.3508079731285591e+01,  4.5071430640991617e+00,  2.3199394181140507e+00,
                9.8658484197036600e-01,  -3.4398135955756350e+00, -3.4400547966379733e+00,
                1.5612856028097871e+01,  3.6617614577493516e+00,  1.0111501174320880e+00,
                2.0807257779604549e+00,  -4.7941550570419755e+00, 4.6990985216199430e+00,
                1.9876039207812369e+01,  -2.8766088932172384e+00, -3.1817503279289938e+00,
                -3.1659549014656618e+00, -1.9575775704046228e+00, 2.4756431632237841e-01,
                9.1393551997910869e+00,  -2.0877085980195806e+00, 1.1185289472237947e+00,
                -3.6050613934795814e+00, -2.0785535146478185e+00, -1.3363815670630830e+00,
                9.5778255802439194e+00};
            // Right-hand side b
            std::array rhs_s{5.7035192278602720e+00, -6.0065243568328057e+00,
                             2.8468876827223211e-01, 1.8482913772408494e+00,
                             -9.0709917456000451e+00};
            // Expected solution x
            std::array x_expected{2.2416392347316735e-01, -1.8161927088157226e-01,
                                  -5.7090973775457202e-02, 1.8212694490425936e-03,
                                  -1.0537576460766331e+00};
            std::array<int, 5> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 5, 5>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 5>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 5>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 5; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-10);
            }
        }
    }

    TEST_CASE("LU solver: identity and diagonal")
    {
        SUBCASE("Identity 3x3")
        {
            // Matrix A
            std::array m_s{
                1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00};
            // Right-hand side b
            std::array rhs_s{1.0000000000000000e+00, 2.0000000000000000e+00,
                             3.0000000000000000e+00};
            // Expected solution x
            std::array x_expected{1.0000000000000000e+00, 2.0000000000000000e+00,
                                  3.0000000000000000e+00};
            std::array<int, 3> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 3; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }

        SUBCASE("Diagonal 4x4")
        {
            // Matrix A
            std::array m_s{
                2.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 3.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 4.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                5.0000000000000000e+00};
            // Right-hand side b
            std::array rhs_s{2.0000000000000000e+00, 6.0000000000000000e+00,
                             1.2000000000000000e+01, 2.0000000000000000e+01};
            // Expected solution x
            std::array x_expected{1.0000000000000000e+00, 2.0000000000000000e+00,
                                  3.0000000000000000e+00, 4.0000000000000000e+00};
            std::array<int, 4> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 4, 4>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 4>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 4>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 4; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }
    }

    TEST_CASE("LU solver: symmetric positive definite")
    {
        SUBCASE("Symmetric positive definite 3x3")
        {
            // Matrix A
            std::array m_s{
                2.4588399791620290e+00, 7.0919365619106389e-01, 4.2507946826226545e-01,
                7.0919365619106389e-01, 1.7595445552758471e+00, 2.3164588162890817e-01,
                4.2507946826226545e-01, 2.3164588162890817e-01, 1.2136829276467180e+00};
            // Right-hand side b
            std::array rhs_s{1.0000000000000000e+00, 2.0000000000000000e+00,
                             3.0000000000000000e+00};
            // Expected solution x
            std::array x_expected{-2.7565413706301506e-01, 9.3307997761344041e-01,
                                  2.3902707320663010e+00};
            std::array<int, 3> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 3; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }

        SUBCASE("Symmetric positive definite 4x4")
        {
            // Matrix A
            std::array m_s{
                2.6423491600080413e+00, 9.4504101139848318e-01, 2.0131778941146559e+00,
                1.6626426068779672e+00, 9.4504101139848318e-01, 1.2074981056618288e+00,
                1.2442162099601155e+00, 1.1045332769598888e+00, 2.0131778941146559e+00,
                1.2442162099601155e+00, 2.0492862612307410e+00, 1.5375459906200610e+00,
                1.6626426068779672e+00, 1.1045332769598888e+00, 1.5375459906200610e+00,
                2.0979381334626752e+00};
            // Right-hand side b
            std::array rhs_s{1.0000000000000000e+00, 1.0000000000000000e+00,
                             1.0000000000000000e+00, 1.0000000000000000e+00};
            // Expected solution x
            std::array x_expected{5.4592409300252620e-01, 1.2206796988264754e+00,
                                  -7.5599242455463356e-01, -4.4608028332488482e-02};
            std::array<int, 4> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 4, 4>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 4>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 4>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 4; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }
    }

    TEST_CASE("LU solver: known solutions")
    {
        SUBCASE("Known solution [1, 2, 3]")
        {
            // Matrix A
            std::array m_s{
                9.5877383705702144e+00,  -2.7897765828645458e+00, -1.3091874747062171e+00,
                -2.0922429623358232e+00, 1.2501928016909789e+01,  2.8204643296019594e+00,
                3.1109436596485125e+00,  -9.9209459226326047e-01, 1.3843046225796916e+01};
            // Right-hand side b
            std::array rhs_s{8.0622780722471532e-02, 3.1373006060289633e+01,
                             4.2655893152512739e+01};
            // Expected solution x
            std::array x_expected{1.0000000000000000e+00, 2.0000000000000000e+00,
                                  3.0000000000000000e+00};
            std::array<int, 3> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 3; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }

        SUBCASE("Random 5x5 diagonally dominant system")
        {
            // Condition number: 1.9605e+00 (well-conditioned)
            // Matrix A
            std::array m_s{
                1.7056125092117890e+01,  -3.8966660991609335e+00, 4.2949510097485781e+00,
                -4.5738788764629373e+00, -2.1614579903147376e+00, 2.3967479659639629e+00,
                1.6935369867040906e+01,  1.5900664645824373e+00,  2.1926699084793092e+00,
                2.0354025260945301e+00,  3.9653220264287877e+00,  -1.4059999751856695e+00,
                1.5899398959869832e+01,  -4.5109850313747949e+00, 4.1885039166734099e+00,
                -4.2177621360009821e+00, -4.2890949082942095e+00, 4.3923776923743134e+00,
                1.9862638734949552e+01,  -2.0091125046113663e+00, 2.9134607902639806e+00,
                -1.2869421104327754e+00, -1.8639960612394334e+00, 3.6043867064799364e+00,
                1.4569326042175085e+01};
            // Right-hand side b
            std::array rhs_s{-4.7098551048808094e+00, 4.6889203568069828e+01,
                             2.8695950988932938e+01, 3.0253318353827959e+01,
                             4.4172758988639825e+01};
            // Expected solution x (verified with NumPy)
            std::array x_expected{6.3872298535535277e-01, 1.9300258496755933e+00,
                                  1.6463940467433373e+00, 1.9938576151440683e+00,
                                  2.7920258063497663e+00};
            std::array<int, 5> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 5, 5>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 5>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 5>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 5; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
            }
        }
    }

    TEST_CASE("LU solver: stiff systems (high condition number)")
    {
        SUBCASE("Diagonal stiff 3x3 (eigenvalue ratio = 1e7)")
        {
            // Condition number: 1.000000e+07
            // Eigenvalue ratio: 1.000000e+07 (min=1e-05, max=100.0)
            // This is a simple diagonal system with known solution [1, 1, 1]
            // Matrix A
            std::array m_s{
                1.0000000000000000e-05, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+02};
            // Right-hand side b
            std::array rhs_s{1.0000000000000000e-05, 1.0000000000000000e+00,
                             1.0000000000000000e+02};
            // Expected solution x
            std::array x_expected{1.0000000000000000e+00, 1.0000000000000000e+00,
                                  1.0000000000000000e+00};
            std::array<int, 3> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution - diagonal system should be very accurate
            for (size_t i = 0; i < 3; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-10);
            }
        }

        SUBCASE("Moderately stiff 4x4 (condition ~ 1e4)")
        {
            // Diagonally scaled system with wide range
            // Matrix A (diagonal matrix with scaling 0.0001, 1, 10, 100)
            std::array m_s{
                1.0000000000000000e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 1.0000000000000000e+01, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                1.0000000000000000e+02};
            // Right-hand side b chosen so solution is [2, 3, 4, 5]
            std::array rhs_s{2.0000000000000000e-04, 3.0000000000000000e+00,
                             4.0000000000000000e+01, 5.0000000000000000e+02};
            // Expected solution x
            std::array x_expected{2.0000000000000000e+00, 3.0000000000000000e+00,
                                  4.0000000000000000e+00, 5.0000000000000000e+00};
            std::array<int, 4> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 4, 4>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 4>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 4>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution
            for (size_t i = 0; i < 4; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-10);
            }
        }

        SUBCASE("Extremely stiff diagonal 5x5 (condition ~ 1e10)")
        {
            // This tests the limits of the solver's numerical stability
            // Eigenvalue ratio 1e-10 to 1.0 = 1e10
            // Matrix A
            std::array m_s{
                1.0000000000000000e-10, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                1.0000000000000000e-05, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                1.0000000000000000e-01, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                1.0000000000000000e+00};
            // Right-hand side b chosen so solution is [1, 1, 1, 1, 1]
            std::array rhs_s{1.0000000000000000e-10, 1.0000000000000000e-05,
                             1.0000000000000000e-01, 1.0000000000000000e+00,
                             1.0000000000000000e+00};
            // Expected solution x
            std::array x_expected{1.0000000000000000e+00, 1.0000000000000000e+00,
                                  1.0000000000000000e+00, 1.0000000000000000e+00,
                                  1.0000000000000000e+00};
            std::array<int, 5> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 5, 5>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 5>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 5>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution with relaxed tolerance for extreme stiffness
            for (size_t i = 0; i < 5; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-8);
            }
        }

        SUBCASE("Large stiff 32x32 diagonal system (condition ~ 1e8)")
        {
            // Diagonal matrix with eigenvalues from 1e-8 to 1.0
            // This tests numerical stability with large, stiff systems
            // Expected solution: x = [1, 2, 3, ..., 32]
            constexpr size_t n = 32;

            // Matrix A (diagonal entries only - rest are zeros)
            std::array<double, 1024> m_s{}; // Initialize all to zero

            // Set diagonal elements
            std::array<double, 32> diagonal{
                1.0000000000000000e-08, 1.8116091942004132e-08, 3.2819278725114711e-08,
                5.9455707085443948e-08, 1.0771050560367691e-07, 1.9512934226359623e-07,
                3.5349811050301095e-07, 6.4040042711972832e-07, 1.1601553017399715e-06,
                2.1017480113324870e-06, 3.8075460212223681e-06, 6.8977853793876580e-06,
                1.2496091412919867e-05, 2.2638034095214464e-05, 4.1011270705513043e-05,
                7.4296395075949503e-05, 1.3459603241553644e-04, 2.4383540982688266e-04,
                4.4173447031400643e-04, 8.0025022781610524e-04, 1.4497406703726315e-03,
                2.6263635276533299e-03, 4.7579443140094140e-03, 8.6195356647530332e-03,
                1.5615230060004965e-02, 2.8288694346259666e-02, 5.1248058769609257e-02,
                9.2841454451947442e-02, 1.6819243248808688e-01, 3.0469895709035050e-01,
                5.5199543212815727e-01, 1.0000000000000000e+00};

            for (size_t i = 0; i < n; ++i) {
                m_s[i * n + i] = diagonal[i];
            }

            // Right-hand side b
            std::array<double, 32> rhs_s{
                1.0000000000000000e-08, 3.6232183884008263e-08, 9.8457836175344139e-08,
                2.3782282834177579e-07, 5.3855252801838451e-07, 1.1707760535815774e-06,
                2.4744867735210768e-06, 5.1232034169578266e-06, 1.0441397715659743e-05,
                2.1017480113324869e-05, 4.1883006233446050e-05, 8.2773424552651893e-05,
                1.6244918836795828e-04, 3.1693247733300252e-04, 6.1516906058269570e-04,
                1.1887423212151920e-03, 2.2881325510641193e-03, 4.3890373768838880e-03,
                8.3929549359661217e-03, 1.6005004556322106e-02, 3.0444554077825262e-02,
                5.7779997608373261e-02, 1.0943271922221652e-01, 2.0686885595407278e-01,
                3.9038075150012413e-01, 7.3550605300275129e-01, 1.3836975867794499e+00,
                2.5995607246545283e+00, 4.8775805421545195e+00, 9.1409687127105155e+00,
                1.7111858395972874e+01, 3.2000000000000000e+01};

            // Expected solution x = [1, 2, 3, ..., 32]
            std::array<double, 32> x_expected{
                1.0000000000000000e+00, 2.0000000000000000e+00, 3.0000000000000000e+00,
                4.0000000000000000e+00, 5.0000000000000000e+00, 6.0000000000000000e+00,
                7.0000000000000000e+00, 8.0000000000000000e+00, 9.0000000000000000e+00,
                1.0000000000000000e+01, 1.1000000000000000e+01, 1.2000000000000000e+01,
                1.3000000000000000e+01, 1.4000000000000000e+01, 1.5000000000000000e+01,
                1.6000000000000000e+01, 1.7000000000000000e+01, 1.8000000000000000e+01,
                1.9000000000000000e+01, 2.0000000000000000e+01, 2.1000000000000000e+01,
                2.2000000000000000e+01, 2.3000000000000000e+01, 2.4000000000000000e+01,
                2.5000000000000000e+01, 2.6000000000000000e+01, 2.7000000000000000e+01,
                2.8000000000000000e+01, 2.9000000000000000e+01, 3.0000000000000000e+01,
                3.1000000000000000e+01, 3.2000000000000000e+01};

            std::array<int, 32> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 32, 32>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 32>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 32>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution - should maintain good accuracy even with high condition
            // number
            for (size_t i = 0; i < 32; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-6);
            }
        }
    }

    TEST_CASE("LU solver: numerically challenging systems")
    {
        SUBCASE("Nearly singular 3x3")
        {
            // Matrix A
            std::array m_s{
                1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00, 1.0000000001000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000001000000e+00};
            // Right-hand side b
            std::array rhs_s{3.0000000000000000e+00, 3.0000000001000000e+00,
                             3.0000000001000000e+00};
            // Expected solution x
            std::array x_expected{1.0000000001000005e+00, 9.9999999900003289e-01,
                                  1.0000000000999670e+00};
            std::array<int, 3> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution with relaxed tolerance
            for (size_t i = 0; i < 3; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-6);
            }
        }

        SUBCASE("Mixed scales 3x3")
        {
            // Matrix A - moderate scaling differences
            // Condition number ~ 1620
            std::array m_s{
                1.0000000000000000e+03, 1.0000000000000000e+00, 1.0000000000000000e+00,
                1.0000000000000000e+00, 1.0000000000000000e-03, 1.0000000000000000e+00,
                1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00};
            // Right-hand side b chosen so solution is [1, 1, 1]
            std::array rhs_s{1.0020000000000000e+03, 2.0010000000000000e+00,
                             3.0000000000000000e+00};
            // Expected solution x (verified with NumPy)
            std::array x_expected{1.0000000000000000e+00, 1.0000000000000002e+00,
                                  9.9999999999999989e-01};
            std::array<int, 3> m_perm_s;

            // Setup mdspan views
            auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
            auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
            auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

            // Solve system
            hd::lu_decomp(m, m_perm);
            hd::lu_backsubs(m, m_perm, rhs);

            // Verify solution with relaxed tolerance
            for (size_t i = 0; i < 3; ++i) {
                CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-10);
            }
        }
    }
}
