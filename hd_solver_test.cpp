#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "fmt/format.h"
#include "fmt/ranges.h"

#include "hd_solver.hpp"

double const eps{1.e-15};

TEST_SUITE("solver:")
{
    TEST_CASE("solver: simple system")
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
}