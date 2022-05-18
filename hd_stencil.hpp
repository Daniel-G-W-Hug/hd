#ifndef HD_STENCIL_H
#define HD_STENCIL_H

#include "hd/hd_functions.hpp" // hd::fact()
#include "hd/hd_solver.hpp"    // hd::lu_decomp(), hd::lu_backsubs()

#include "mdspan/mdspan.hpp"

#include "fmt/format.h"
#include "fmt/ranges.h"

#include <cmath>     // std::pow(), std::abs()
#include <stdexcept> //std::invalid_argument
#include <vector>

namespace hd // Namespace HD to define my types for numerical computation
{

enum class stencil_lhs
{
    f1, // f'  terms considered to be on lhs of finite difference
    f2  // f'' terms considered to be on lhs of finite difference
};

struct stencil_t
{
    // after calling the ctor, the "output values" can be used
    stencil_t(double x0, stencil_lhs lhs_t, std::vector<double> xf0, std::vector<double> xf1, std::vector<double> xf2) :
        x0{x0}, lhs_t{lhs_t}, xf0{xf0}, xf1{xf1}, xf2{xf2}
    {

        // consistency checks
        if ((nf1() == 0 && nf2() == 0) || n() < 3 ||
            (nf1() == 0 && lhs_t == stencil_lhs::f1) ||
            (nf2() == 0 && lhs_t == stencil_lhs::f2))
        {
            throw std::invalid_argument("Inconsistent specification of stencil in ctor of hd::stencil_t.");
        }

        // reserve space for weights (does not change wfx.size())
        wf0.reserve(xf0.size());
        wf1.reserve(xf1.size());
        wf2.reserve(xf2.size());

        // fmt::print("wf0.size()={}, wf0.capacity()={}\n", wf0.size(), wf0.capacity());
        // fmt::print("wf1.size()={}, wf1.capacity()={}\n", wf1.size(), wf1.capacity());
        // fmt::print("wf2.size()={}, wf2.capacity()={}\n", wf2.size(), wf2.capacity());

        // provide the weights, order and truncation error for further use
        calc_stencil();
    }

    // input values provided by ctor
    const double x0;         // development point of stencil
                             // (should be within or at least close to coordinates of points)
    const stencil_lhs lhs_t; // either f1 or f2 terms on lhs, all other terms considered to be on rhs

    const std::vector<double> xf0; // coordinates of points for f
    const std::vector<double> xf1; // coordinates of points for f'
    const std::vector<double> xf2; // coordinates of points for f''

    // output values filled in after call of calc_stencil()
    std::vector<double> wf0; // weights of points for f
    std::vector<double> wf1; // weights of points for f'
    std::vector<double> wf2; // weights of points for f''

    int order;        // order of fd stencil
    double trunc_err; // truncation error as factor in front of highest neglected term

    // helpers for number of points
    int nf0() const { return xf0.size(); }          // number of points for f
    int nf1() const { return xf1.size(); }          // number of points for f'
    int nf2() const { return xf2.size(); }          // number of points for f''
    int n() const { return nf0() + nf1() + nf2(); } // total number of points

  private:

    void calc_stencil();
};

void stencil_t::calc_stencil()
{

    // reserve memory for matrix, permutation and rhs vector and initialize with 0.0
    std::vector<double> mem_matrix(n() * n(), 0.0);
    std::vector<int> mem_perm(n(), 0.0);
    std::vector<double> mem_rhs(n(), 0.0);

    // create views
    mdspan matrix{mem_matrix.data(), n(), n()};
    mdspan perm{mem_perm.data(), n()};
    mdspan rhs{mem_rhs.data(), n()};

    // setup column indices (i.e. begin/end indices for f, f', f'')
    int j0b, j0e, j1b, j1e, j2b, j2e, col;

    // f
    j0b = 0;
    j0e = nf0() - 1;
    col = j0e;

    // f'
    if (nf1() > 0)
    {
        j1b = col + 1;
        j1e = col + nf1();
        col += nf1();
    }

    // f''
    if (nf2() > 0)
    {
        j2b = col + 1;
        j2e = col + nf2();
        // col += nf2(); // just needed for further extension with f'''
    }

    // fmt::print("nf0()={}, nf1()={}, nf2()={}, n()={}\n", nf0(), nf1(), nf2(), n());
    // fmt::print("j0b={}, j0e={}, j1b={}, j1e={}, j2b={}, j2e={}\n", j0b, j0e, j1b, j1e, j2b, j2e);

    // setup standard matrix columnwise
    // unwanted terms of series expansion of lhs are moved to rhs (sfact)

    // f
    for (int j = j0b; j <= j0e; ++j)
    {
        matrix(0, j) = 1.0;
        for (int i = 1; i < n(); ++i)
            matrix(i, j) = std::pow(xf0[j - j0b] - x0, i) / hd::fact(i);
    }

    double sfact; // used to locate terms on lhs (-1.0) or on rhs (+1.0)

    // f'
    if (nf1() > 0)
    {
        if (lhs_t == stencil_lhs::f1)
        {
            // put terms on lhs
            sfact = -1.0;
        }
        else
        {
            // put terms on rhs
            sfact = 1.0;
        }
        for (int j = j1b; j <= j1e; ++j)
        {
            matrix(0, j) = 0.0;
            matrix(1, j) = 1.0;
            for (int i = 2; i < n(); ++i)
                matrix(i, j) = sfact * std::pow(xf1[j - j1b] - x0, i - 1) / hd::fact(i - 1);
        }
    }

    // f''
    if (nf2() > 0)
    {
        if (lhs_t == stencil_lhs::f2)
        {
            // put terms on lhs
            sfact = -1.0;
        }
        else
        {
            // put terms on rhs
            sfact = 1.0;
        }
        for (int j = j2b; j <= j2e; ++j)
        {
            matrix(0, j) = 0.0;
            matrix(1, j) = 0.0;
            matrix(2, j) = 1.0;
            for (int i = 3; i < n(); ++i)
                matrix(i, j) = sfact * std::pow(xf2[j - j2b] - x0, i - 2) / hd::fact(i - 2);
        }
    }

    // setup rhs
    if (lhs_t == stencil_lhs::f1)
    {
        rhs(1) = 1.0;
    }

    if (lhs_t == stencil_lhs::f2)
    {
        rhs(2) = 1.0;
    }

    // normalization: replace last equation with normalization condition (sum of coefficients on lhs = 1)
    //                i.e. set coefficients of primary derivative to 1.0 in the last equation (normalization)
    //                and set them to 0.0 in the corresponding matrix line
    //                (remove them from the rhs of the standard system)
    for (int j = 0; j < n(); ++j)
        matrix(n() - 1, j) = 0.0;

    rhs(n() - 1) = 1.0;

    if (lhs_t == stencil_lhs::f1)
    {
        for (int j = j1b; j <= j1e; ++j)
        {
            matrix(n() - 1, j) = 1.0;
            matrix(1, j) = 0.0;
        }
    }
    if (lhs_t == stencil_lhs::f2)
    {
        for (int j = j2b; j <= j2e; ++j)
        {
            matrix(n() - 1, j) = 1.0;
            matrix(2, j) = 0.0;
        }
    }

    // solve system
    hd::lu_decomp(matrix, perm);
    hd::lu_backsubs(matrix, perm, rhs); // weights are now on rhs

    // assign weights to output vectors
    // f
    for (int j = j0b; j <= j0e; ++j)
        wf0.push_back(rhs(j));
    // f'
    if (nf1() > 0)
    {
        for (int j = j1b; j <= j1e; ++j)
            wf1.push_back(rhs(j));
    }
    // f''
    if (nf2() > 0)
    {
        for (int j = j2b; j <= j2e; ++j)
            wf2.push_back(rhs(j));
    }

    // compute order and truncation error

    // f
    for (int i = nf0(); i <= n(); ++i)
    {
        double sumte = 0.0;
        for (int j = j0b; j <= j0e; ++j)
            sumte += std::pow(xf0[j - j0b] - x0, i) / hd::fact(i) * rhs(j);

        // f'
        if (nf1() > 0)
        {
            if (lhs_t == stencil_lhs::f1)
            {
                // put terms on lhs
                sfact = -1.0;
            }
            else
            {
                // put terms on rhs
                sfact = 1.0;
            }
            for (int j = j1b; j <= j1e; ++j)
                sumte += std::pow(xf1[j - j1b] - x0, i - 1) / hd::fact(i - 1) * rhs(j);
        }

        // f''
        if (nf2() > 0)
        {
            if (lhs_t == stencil_lhs::f2)
            {
                // put terms on lhs
                sfact = -1.0;
            }
            else
            {
                // put terms on rhs
                sfact = 1.0;
            }
            for (int j = j2b; j <= j2e; ++j)
                sumte += std::pow(xf2[j - j2b] - x0, i - 2) / hd::fact(i - 2) * rhs(j);
        }

        double eps = 1.0e-6;
        if (std::abs(sumte) > eps)
        {
            trunc_err = sumte;

            if (lhs_t == stencil_lhs::f1)
            {
                order = i - 1;
            }
            if (lhs_t == stencil_lhs::f2)
            {
                order = i - 2;
            }
        }
    }
}

} // namespace hd

#endif // HD_STENCIL_H