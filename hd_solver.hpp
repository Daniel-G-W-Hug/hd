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

// use branch "single-header" from mdspan github
//
// To try using subscript operator comment in macro below
// the header will by default also check for the feature macro, and enable it
// defining the macro to 0 will overwrite the automatic setting
// x86-64 clang (experimental auto NSDMI) supports the operator, but you need
// to explicitly comment in below macro
// #define MDSPAN_USE_BRACKET_OPERATOR 1

// To force enable operator() comment in the macro below
// You can enable both at the same time.
// #define MDSPAN_USE_PAREN_OPERATOR 1

#include "mdspan/mdspan.hpp"

#include <cmath>
#include <iostream>
#include <vector>

// make mdspan less verbose
using namespace Kokkos;

namespace hd { // Namespace hd to define my types for numerical computation

void lu_decomp(mdspan<double, dextents<size_t, 2>> a,
               mdspan<int, dextents<size_t, 1>> perm);
void lu_backsubs(mdspan<double const, dextents<size_t, 2>> a,
                 mdspan<int const, dextents<size_t, 1>> perm,
                 mdspan<double, dextents<size_t, 1>> b);

//-----------------------------------------------------------------------------
// Solver error handling
//-----------------------------------------------------------------------------
struct Solver_error {
    std::string name;
    Solver_error(const char* q) : name(q) {}
    Solver_error(std::string n) : name(n) {}
};

//-----------------------------------------------------------------------------

inline void solver_error_msg(const char* p) { throw Solver_error(p); }

void lu_decomp(mdspan<double, dextents<size_t, 2>> a,
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
            if (abs(a[i, j]) > aamax) aamax = abs(a[i, j]);
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
            dum = vv[i] * abs(sum);
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

} // ludecomp()

void lu_backsubs(mdspan<double const, dextents<size_t, 2>> a,
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

        solver_error_msg("hd::lu_decomp(): unsymmetric matrix, permututation vector size "
                         "or right hand side size incompatible.");
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

} // lubacksubs()

} // namespace hd

#endif // HD_SOLVER_H