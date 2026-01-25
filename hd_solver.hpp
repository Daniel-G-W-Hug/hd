#ifndef HD_SOLVER_H
#define HD_SOLVER_H

// implementation of solvers for small systems a*x = b (gaussian elimination by LU
// factorization of a)
//
// Usage:
//
// 1.) LU decomposition of matrix
//
// hd::lu_decomp(m, m_perm);
//
// 2.) solution by backsubstition of rhs => solution is returned in rhs
//     (can be repeated with many different rhs vectors for the same matrix)
//
// hd::lu_backsubs(m, m_perm, rhs);
//
// mdspan selection:
// - By default, uses C++23 standard <mdspan> if available (__cpp_lib_mdspan)
// - Define HD_SOLVER_USE_KOKKOS_MDSPAN to force use of Kokkos mdspan implementation
// - Define HD_SOLVER_USE_STD_MDSPAN to force use of C++23 standard <mdspan>

// Include <version> for feature test macros (C++20+)
#if __has_include(<version>)
    #include <version>
#endif

// mdspan selection logic
#if defined(HD_SOLVER_USE_KOKKOS_MDSPAN)
    // User explicitly requested Kokkos mdspan
    #include "mdspan/mdspan.hpp"
    #define HD_SOLVER_MDSPAN_NS Kokkos
#elif defined(HD_SOLVER_USE_STD_MDSPAN) || (defined(__cpp_lib_mdspan) && __cpp_lib_mdspan >= 202207L)
    // Use C++23 standard mdspan (explicitly requested or auto-detected)
    #include <mdspan>
    #define HD_SOLVER_MDSPAN_NS std
#else
    // Fallback to Kokkos mdspan for older compilers
    #include "mdspan/mdspan.hpp"
    #define HD_SOLVER_MDSPAN_NS Kokkos
#endif

#include <cmath>
#include <stdexcept>
#include <vector>

namespace hd { // Namespace hd to define my types for numerical computation

// Use selected mdspan namespace
using HD_SOLVER_MDSPAN_NS::dextents;
using HD_SOLVER_MDSPAN_NS::extents;
using HD_SOLVER_MDSPAN_NS::mdspan;

//-----------------------------------------------------------------------------
// Solver error handling
//-----------------------------------------------------------------------------
struct Solver_error {
    std::string name;
    Solver_error(char const* q) : name(q) {}
    Solver_error(std::string n) : name(n) {}
};

//-----------------------------------------------------------------------------

inline void solver_error_msg(char const* p) { throw Solver_error(p); }

//-----------------------------------------------------------------------------
// LU decomposition
//-----------------------------------------------------------------------------

inline void lu_decomp(mdspan<double, dextents<size_t, 2>> a,
                      mdspan<int, dextents<size_t, 1>> perm)
{
    /* LU decomposition of matrix a (handed back on a)
       perm is the permutation vector in case of line exchange (pivot elements)
    */

    // check fitness of matrix and permutation vector
    if (a.extent(0) != a.extent(1) || a.extent(0) != perm.extent(0)) {

        solver_error_msg("hd::lu_decomp(): unsymmetric matrix or permututation vector "
                         "size incompatible.");
    };

    constexpr double TINY = 1.e-20;
    int ubound = static_cast<int>(a.extent(0)) - 1; // highest valid index (=upper boundary)

    // helper for scaling the matrix rows
    std::vector<double> vv(a.extent(0));

    // fill in scaling vector
    for (int i = 0; i <= ubound; ++i) {
        double aamax = 0.;
        for (int j = 0; j <= ubound; ++j) {
            if (std::abs(a[i, j]) > aamax) aamax = std::abs(a[i, j]);
        }
        if (aamax == 0.) solver_error_msg("hd::lu_decomp(): singular matrix.");
        vv[i] = 1. / aamax;
    }

    // LU decomposition
    double sum, aamax, dum;
    int imax = 0; // Initialize to avoid MSVC warning
    for (int j = 0; j <= ubound; ++j) {
        if (j > 0) {
            for (int i = 0; i <= j - 1; ++i) {
                sum = a[i, j];
                if (i > 0) {
                    for (int k = 0; k <= i - 1; ++k)
                        sum -= a[i, k] * a[k, j];
                    a[i, j] = sum;
                }
            }
        }
        aamax = 0.;
        for (int i = j; i <= ubound; ++i) {
            sum = a[i, j];
            if (j > 0) {
                for (int k = 0; k <= j - 1; ++k)
                    sum -= a[i, k] * a[k, j];
                a[i, j] = sum;
            }
            dum = vv[i] * std::abs(sum);
            if (dum >= aamax) {
                imax = i;
                aamax = dum;
            }
        }
        if (j != imax) {
            for (int k = 0; k <= ubound; ++k) {
                dum = a[imax, k];
                a[imax, k] = a[j, k];
                a[j, k] = dum;
            }
            vv[imax] = vv[j];
        }
        perm[j] = imax;
        if (j != ubound) {
            if (a[j, j] == 0.) a[j, j] = TINY;
            dum = 1. / a[j, j];
            for (int i = j + 1; i <= ubound; ++i)
                a[i, j] *= dum;
        }
    }
    if (a[ubound, ubound] == 0.) a[ubound, ubound] = TINY;

} // lu_decomp()

//-----------------------------------------------------------------------------
// Backward substitution
//-----------------------------------------------------------------------------

inline void lu_backsubs(mdspan<double const, dextents<size_t, 2>> a,
                        mdspan<int const, dextents<size_t, 1>> perm,
                        mdspan<double, dextents<size_t, 1>> b)
{
    /*
    backward substitution: a is the LU-decomposed matrix as provided by lu_decomp()
    perm is the corresponding permutation vector provided by lu_decomp()

    lu_decomp() must be called once before lu_backsubs

    b is the right hand side of the equation a*x = b

    The solution vector x will be returned on b

    lu_backsubs() can be used for arbitrarily many different right hand side vectors
    */

    // check fitness of matrix, permutation vector and right hand side
    if (a.extent(0) != a.extent(1) || a.extent(0) != perm.extent(0) ||
        a.extent(0) != b.extent(0)) {

        solver_error_msg("hd::lu_backsubs(): unsymmetric matrix, permututation vector "
                         "size or right hand side size incompatible.");
    };

    int ubound = static_cast<int>(a.extent(0)) - 1; // highest valid index (=upper boundary)

    double sum;
    int ll;

    int ii = -1; // never occurring index as indicator for first loop
    for (int i = 0; i <= ubound; ++i) {
        ll = perm[i];
        sum = b[ll];
        b[ll] = b[i];
        if (ii != -1) {
            for (int j = ii; j <= i - 1; ++j)
                sum -= a[i, j] * b[j];
        }
        else if (sum != 0.) ii = i;
        b[i] = sum;
    }

    for (int i = ubound; i >= 0; --i) {
        sum = b[i];
        if (i < ubound) {
            for (int j = i + 1; j <= ubound; ++j)
                sum -= a[i, j] * b[j];
        }
        b[i] = sum / a[i, i];
    }

} // lu_backsubs()

} // namespace hd

#endif // HD_SOLVER_H
