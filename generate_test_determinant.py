#!/usr/bin/env python3
"""
Generate test matrices for determinant testing.

This script generates random matrices with known determinants
to be used as test cases in hd_determinant_test.cpp.

Usage:
    python3 generate_test_matrices.py
"""

import numpy as np


def print_matrix(A, name, det_value):
    """Print a matrix in C++ vector format."""
    n = A.shape[0]
    print(f"\n// {name}")
    print(f"// Determinant: {det_value}")
    print("std::vector<std::vector<double>> A = {")
    for i, row in enumerate(A):
        comma = "," if i < n - 1 else ""
        print("    {" + ", ".join(f"{x:.6f}" for x in row) + "}" + comma)
    print("};")
    print(f"CHECK(det(A) == doctest::Approx({det_value:.10f}).epsilon(1e-6));")


def print_mdspan_matrix(A, name, det_value):
    """Print a matrix in C++ mdspan format."""
    n = A.shape[0]
    print(f"\n// {name}")
    print(f"// Determinant: {det_value}")
    print("double data[] = {")
    for i in range(n):
        for j in range(n):
            comma = "," if not (i == n-1 and j == n-1) else ""
            if j == 0:
                print("    ", end="")
            print(f"{A[i,j]:.6f}{comma}", end="")
            if j == n-1:
                print()
            else:
                print(" ", end="")
    print("};")
    print(f"Kokkos::mdspan<double, Kokkos::dextents<size_t, 2>> A(data, {n}, {n});")
    print(f"CHECK(det<double>(A) == doctest::Approx({det_value:.10f}).epsilon(1e-6));")


def generate_test_matrices():
    """Generate various test matrices."""

    print("="*80)
    print("Test Matrices for hd_determinant_test.cpp")
    print("="*80)

    # 2x2 matrices
    print("\n" + "="*80)
    print("2x2 Matrices")
    print("="*80)

    A = np.array([[1, 2], [3, 4]])
    det_A = np.linalg.det(A)
    print_matrix(A, "Simple 2x2 matrix", det_A)

    A = np.array([[2, 3], [1, 4]])
    det_A = np.linalg.det(A)
    print_matrix(A, "Integer 2x2 matrix", det_A)

    # 3x3 matrices
    print("\n" + "="*80)
    print("3x3 Matrices")
    print("="*80)

    A = np.array([[1, 2, 3], [0, 1, 4], [5, 6, 0]])
    det_A = np.linalg.det(A)
    print_matrix(A, "Simple 3x3 matrix", det_A)

    A = np.array([[0, 2, 6], [1, 8, 4], [5, 2, 7]])
    det_A = np.linalg.det(A)
    print_matrix(A, "User example 3x3 matrix", det_A)

    # 4x4 matrices
    print("\n" + "="*80)
    print("4x4 Matrices")
    print("="*80)

    A = np.array([[1, 2, 0, 1],
                  [3, 1, 2, 0],
                  [0, 1, 1, 2],
                  [2, 0, 1, 1]])
    det_A = np.linalg.det(A)
    print_matrix(A, "General 4x4 matrix", det_A)

    A = np.diag([2, 3, 4, 5])
    det_A = np.linalg.det(A)
    print_matrix(A, "Diagonal 4x4 matrix", det_A)

    # 5x5 matrices
    print("\n" + "="*80)
    print("5x5 Matrices")
    print("="*80)

    np.random.seed(42)  # For reproducibility
    A = np.random.uniform(-1, 1, (5, 5))
    det_A = np.linalg.det(A)
    print_matrix(A, "Random 5x5 matrix with entries in [-1, 1]", det_A)

    # Another random 5x5
    np.random.seed(123)
    A = np.random.uniform(-1, 1, (5, 5))
    det_A = np.linalg.det(A)
    print_matrix(A, "Another random 5x5 matrix", det_A)

    # 6x6 matrix
    print("\n" + "="*80)
    print("6x6 Matrices")
    print("="*80)

    np.random.seed(999)
    A = np.random.uniform(-1, 1, (6, 6))
    det_A = np.linalg.det(A)
    print_matrix(A, "Random 6x6 matrix with entries in [-1, 1]", det_A)

    # Special matrices
    print("\n" + "="*80)
    print("Special Matrices")
    print("="*80)

    # Upper triangular
    A = np.array([[1, 2, 3, 4],
                  [0, 1, 2, 3],
                  [0, 0, 1, 2],
                  [0, 0, 0, 1]])
    det_A = np.linalg.det(A)
    print_matrix(A, "Upper triangular 4x4 (det = product of diagonal)", det_A)

    # Orthogonal matrix (det = ±1)
    theta = np.pi / 4
    A = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    det_A = np.linalg.det(A)
    print_matrix(A, "Rotation matrix (det = 1)", det_A)

    print("\n" + "="*80)


def generate_singular_matrices():
    """Generate singular matrices (determinant = 0)."""

    print("\n" + "="*80)
    print("Singular Matrices (det = 0)")
    print("="*80)

    # Row of zeros
    A = np.array([[1, 2, 3],
                  [0, 0, 0],
                  [4, 5, 6]])
    det_A = np.linalg.det(A)
    print_matrix(A, "3x3 with row of zeros", det_A)

    # Linearly dependent rows
    A = np.array([[1, 2, 3, 4],
                  [2, 4, 6, 8],
                  [1, 1, 1, 1],
                  [0, 1, 2, 3]])
    det_A = np.linalg.det(A)
    print_matrix(A, "4x4 with linearly dependent rows", det_A)


if __name__ == "__main__":
    generate_test_matrices()
    generate_singular_matrices()

    print("\n" + "="*80)
    print("Note: Copy relevant matrices to hd_determinant_test.cpp")
    print("="*80)
