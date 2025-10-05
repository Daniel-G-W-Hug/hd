#ifndef HD_DETERMINANT_H
#define HD_DETERMINANT_H

#include "hd_solver.hpp"

#include <stdexcept>
#include <vector>

// make mdspan less verbose
using namespace Kokkos;

namespace hd {

// Determinant using vector of vectors with LU decomposition
template <typename T> T det(const std::vector<std::vector<T>>& A)
{
    const size_t n = A.size();

    if (n == 0 || A[0].size() != n) {
        throw std::invalid_argument("det: Matrix must be square and non-empty");
    }

    // Copy matrix to mdspan-compatible storage
    std::vector<double> data(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            data[i * n + j] = static_cast<double>(A[i][j]);
        }
    }

    mdspan<double, dextents<size_t, 2>> a(data.data(), n, n);
    std::vector<int> perm_data(n);
    mdspan<int, dextents<size_t, 1>> perm(perm_data.data(), n);

    try {
        lu_decomp(a, perm);
    }
    catch (const Solver_error&) {
        return T(0); // Singular matrix has determinant 0
    }

    // Determinant = product of diagonal * sign from permutations
    double result = 1.0;
    for (size_t i = 0; i < n; ++i) {
        result *= a[i, i];
    }

    // Count row swaps to determine sign
    int swaps = 0;
    for (size_t i = 0; i < n; ++i) {
        if (perm[i] != static_cast<int>(i)) {
            ++swaps;
        }
    }

    return static_cast<T>((swaps % 2 == 0) ? result : -result);
}

// Determinant using mdspan with LU decomposition
template <typename T, typename Extents, typename LayoutPolicy, typename AccessorPolicy>
T det(mdspan<T, Extents, LayoutPolicy, AccessorPolicy> A)
{
    const size_t n = A.extent(0);

    if (n != A.extent(1)) {
        throw std::invalid_argument("det: Matrix must be square");
    }

    // Copy to double storage for lu_decomp
    std::vector<double> data(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            data[i * n + j] = static_cast<double>(A[i, j]);
        }
    }

    mdspan<double, dextents<size_t, 2>> a(data.data(), n, n);
    std::vector<int> perm_data(n);
    mdspan<int, dextents<size_t, 1>> perm(perm_data.data(), n);

    try {
        lu_decomp(a, perm);
    }
    catch (const Solver_error&) {
        return T(0); // Singular matrix
    }

    // Product of diagonal elements
    double result = 1.0;
    for (size_t i = 0; i < n; ++i) {
        result *= a[i, i];
    }

    // Sign correction from row swaps
    int swaps = 0;
    for (size_t i = 0; i < n; ++i) {
        if (perm[i] != static_cast<int>(i)) {
            ++swaps;
        }
    }

    return static_cast<T>((swaps % 2 == 0) ? result : -result);
}

} // namespace hd

#endif // HD_DETERMINANT_H
