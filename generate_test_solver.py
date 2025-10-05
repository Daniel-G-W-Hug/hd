#!/usr/bin/env python3
"""
Generate test cases for LU decomposition solver testing.

This script generates random linear systems Ax = b with known solutions
to be used as test cases in hd_solver_test.cpp. It uses NumPy's solver
as a reference implementation.

Usage:
    python3 generate_test_solver.py
"""

import numpy as np


def print_solver_test_case(A, b, x_expected, name, epsilon="1e-12"):
    """Print a solver test case in C++ doctest format with mdspan."""
    n = A.shape[0]

    print(f"\n        SUBCASE(\"{name}\")")
    print("        {")

    # Print matrix A in row-major format
    print("            // Matrix A")
    print("            std::array m_s{", end="")
    elements = []
    for i in range(n):
        for j in range(n):
            elements.append(f"{A[i, j]:.16e}")
    print(", ".join(elements) + "};")

    # Print RHS vector b
    print("            // Right-hand side b")
    print("            std::array rhs_s{", end="")
    elements = [f"{b[i]:.16e}" for i in range(n)]
    print(", ".join(elements) + "};")

    # Print expected solution
    print("            // Expected solution x")
    print("            std::array x_expected{", end="")
    elements = [f"{x_expected[i]:.16e}" for i in range(n)]
    print(", ".join(elements) + "};")

    # Print permutation array
    print(f"            std::array<int, {n}> m_perm_s;")
    print()

    # Setup mdspan views
    print("            // Setup mdspan views")
    print(f"            auto m = mdspan<double, extents<size_t, {n}, {n}>>(m_s.data());")
    print(f"            auto rhs = mdspan<double, extents<size_t, {n}>>(rhs_s.data());")
    print(f"            auto m_perm = mdspan<int, extents<size_t, {n}>>(m_perm_s.data());")
    print()

    # Solve system
    print("            // Solve system")
    print("            hd::lu_decomp(m, m_perm);")
    print("            hd::lu_backsubs(m, m_perm, rhs);")
    print()

    # Check solution
    print("            // Verify solution")
    print(f"            for (size_t i = 0; i < {n}; ++i) {{")
    print(f"                CHECK(std::abs(rhs[i] - x_expected[i]) < {epsilon});")
    print("            }")

    print("        }")


def print_matrix_info(A, name):
    """Print diagnostic information about a matrix."""
    print(f"\n// {name}")
    print(f"// Matrix size: {A.shape[0]}x{A.shape[1]}")
    print(f"// Determinant: {np.linalg.det(A):.6e}")
    cond = np.linalg.cond(A)
    print(f"// Condition number: {cond:.6e}")
    eigenvalues = np.linalg.eigvals(A)
    print(f"// Eigenvalues: min={np.min(np.abs(eigenvalues)):.6e}, max={np.max(np.abs(eigenvalues)):.6e}")
    print(f"// Eigenvalue ratio: {np.max(np.abs(eigenvalues))/np.min(np.abs(eigenvalues)):.6e}")


def generate_basic_tests():
    """Generate basic test cases for various matrix sizes."""

    print("="*80)
    print("Basic Test Cases")
    print("="*80)

    # 2x2 system
    A = np.array([[2.0, 1.0],
                  [1.0, 3.0]])
    b = np.array([5.0, 6.0])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Simple 2x2 system")
    print_solver_test_case(A, b, x, "Simple 2x2 system")

    # 3x3 system
    A = np.array([[1.0, 2.0, 3.0],
                  [0.0, 4.0, 1.0],
                  [0.0, 0.0, 1.0]])
    b = np.array([1.0, 1.0, 1.0])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Upper triangular 3x3 system")
    print_solver_test_case(A, b, x, "Upper triangular 3x3 system")

    # 4x4 system - use a different matrix that's not singular
    A = np.array([[2.0, 1.0, 1.0, 0.0],
                  [4.0, 3.0, 3.0, 1.0],
                  [8.0, 7.0, 9.0, 5.0],
                  [6.0, 7.0, 9.0, 8.0]])
    b = np.array([1.0, 2.0, 3.0, 4.0])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "General 4x4 system")
    print_solver_test_case(A, b, x, "General 4x4 system")

    # 5x5 system
    np.random.seed(42)
    A_temp = np.random.uniform(-5, 5, (5, 5))
    # Make it diagonally dominant for stability
    A = A_temp.copy()
    for i in range(5):
        A[i, i] = np.sum(np.abs(A_temp[i, :])) + 1.0
    b = np.random.uniform(-10, 10, 5)
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Random 5x5 diagonally dominant system")
    print_solver_test_case(A, b, x, "Random 5x5 diagonally dominant system", epsilon="1e-10")


def generate_identity_tests():
    """Generate tests with identity and diagonal matrices."""

    print("\n" + "="*80)
    print("Identity and Diagonal Matrix Tests")
    print("="*80)

    # Identity 3x3
    A = np.eye(3)
    b = np.array([1.0, 2.0, 3.0])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Identity 3x3")
    print_solver_test_case(A, b, x, "Identity 3x3")

    # Diagonal 4x4
    A = np.diag([2.0, 3.0, 4.0, 5.0])
    b = np.array([2.0, 6.0, 12.0, 20.0])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Diagonal 4x4")
    print_solver_test_case(A, b, x, "Diagonal 4x4")


def generate_symmetric_positive_definite_tests():
    """Generate symmetric positive definite systems."""

    print("\n" + "="*80)
    print("Symmetric Positive Definite Systems")
    print("="*80)

    # 3x3 SPD
    # Generate SPD matrix as A = B^T * B
    np.random.seed(100)
    B = np.random.rand(3, 3)
    A = B.T @ B + np.eye(3)  # Add identity to ensure positive definiteness
    b = np.array([1.0, 2.0, 3.0])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Symmetric positive definite 3x3")
    print_solver_test_case(A, b, x, "Symmetric positive definite 3x3")

    # 4x4 SPD
    np.random.seed(200)
    B = np.random.rand(4, 4)
    A = B.T @ B + 0.1 * np.eye(4)
    b = np.ones(4)
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Symmetric positive definite 4x4")
    print_solver_test_case(A, b, x, "Symmetric positive definite 4x4")


def generate_stiff_system_tests():
    """Generate stiff systems with high condition numbers."""

    print("\n" + "="*80)
    print("Stiff System Tests (High Condition Number)")
    print("="*80)

    # Create a 3x3 stiff system with condition number ~ 1e5
    # Use diagonal matrix with large eigenvalue spread
    eigenvalues = np.array([1e-5, 1.0, 1.0])
    # Create orthogonal matrix for eigenvectors
    np.random.seed(500)
    Q, _ = np.linalg.qr(np.random.rand(3, 3))
    # A = Q * Lambda * Q^T
    A = Q @ np.diag(eigenvalues) @ Q.T
    b = np.array([1.0, 2.0, 3.0])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Stiff 3x3 system (condition number ~ 1e5)")
    print_solver_test_case(A, b, x, "Stiff 3x3 system (condition ~ 1e5)", epsilon="1e-8")

    # Create a 4x4 stiff system with condition number ~ 1e6
    eigenvalues = np.array([1e-6, 0.01, 0.1, 1.0])
    np.random.seed(501)
    Q, _ = np.linalg.qr(np.random.rand(4, 4))
    A = Q @ np.diag(eigenvalues) @ Q.T
    b = np.ones(4)
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Stiff 4x4 system (condition number ~ 1e6)")
    print_solver_test_case(A, b, x, "Stiff 4x4 system (condition ~ 1e6)", epsilon="1e-6")

    # Create a 5x5 extremely stiff system with condition number ~ 1e7
    eigenvalues = np.array([1e-7, 0.001, 0.01, 0.1, 1.0])
    np.random.seed(502)
    Q, _ = np.linalg.qr(np.random.rand(5, 5))
    A = Q @ np.diag(eigenvalues) @ Q.T
    b = np.arange(1.0, 6.0)
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Extremely stiff 5x5 system (condition number ~ 1e7)")
    print_solver_test_case(A, b, x, "Extremely stiff 5x5 system (condition ~ 1e7)", epsilon="1e-5")

    # Create a diagonal stiff system (easier to verify)
    A = np.diag([1e-5, 1.0, 100.0])
    b = np.array([1e-5, 1.0, 100.0])  # Solution will be [1, 1, 1]
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Diagonal stiff 3x3 (eigenvalue ratio = 1e7)")
    print_solver_test_case(A, b, x, "Diagonal stiff 3x3 (eigenvalue ratio = 1e7)", epsilon="1e-10")


def generate_tricky_systems():
    """Generate systems that might be numerically challenging."""

    print("\n" + "="*80)
    print("Numerically Challenging Systems")
    print("="*80)

    # Nearly singular (but still invertible) system
    A = np.array([[1.0, 1.0, 1.0],
                  [1.0, 1.0 + 1e-10, 1.0],
                  [1.0, 1.0, 1.0 + 1e-10]])
    b = np.array([3.0, 3.0 + 1e-10, 3.0 + 1e-10])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Nearly singular 3x3")
    print_solver_test_case(A, b, x, "Nearly singular 3x3", epsilon="1e-6")

    # System with mixed scales
    A = np.array([[1e6, 1.0, 1.0],
                  [1.0, 1e-6, 1.0],
                  [1.0, 1.0, 1.0]])
    b = np.array([1e6 + 2.0, 1.0 + 1e-6 + 1.0, 3.0])
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "Mixed scales 3x3")
    print_solver_test_case(A, b, x, "Mixed scales 3x3", epsilon="1e-6")


def generate_known_solution_tests():
    """Generate systems where the solution is known by construction."""

    print("\n" + "="*80)
    print("Known Solution Tests (x = [1, 2, 3, ...])")
    print("="*80)

    # 3x3 with known solution [1, 2, 3]
    np.random.seed(300)
    A_temp = np.random.rand(3, 3) * 10 - 5  # Random matrix in [-5, 5]
    # Make diagonally dominant
    A = A_temp.copy()
    for i in range(3):
        A[i, i] = np.sum(np.abs(A_temp[i, :])) + 5.0
    x_true = np.array([1.0, 2.0, 3.0])
    b = A @ x_true  # Compute RHS from known solution
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "System with known solution [1, 2, 3]")
    print_solver_test_case(A, b, x_true, "Known solution [1, 2, 3]")

    # 5x5 with known solution [1, 2, 3, 4, 5]
    np.random.seed(301)
    A_temp = np.random.rand(5, 5) * 10 - 5
    A = A_temp.copy()
    for i in range(5):
        A[i, i] = np.sum(np.abs(A_temp[i, :])) + 5.0
    x_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    b = A @ x_true
    x = np.linalg.solve(A, b)
    print_matrix_info(A, "System with known solution [1, 2, 3, 4, 5]")
    print_solver_test_case(A, b, x_true, "Known solution [1, 2, 3, 4, 5]")


def print_usage_instructions():
    """Print instructions for using the generated test cases."""
    print("\n" + "="*80)
    print("Usage Instructions")
    print("="*80)
    print("""
To use these test cases:

1. Run this script:
   python3 generate_test_solver.py > new_test_cases.txt

2. Copy the relevant test subcases into hd_solver_test.cpp

3. Organize them into appropriate TEST_CASE sections, such as:
   - TEST_CASE("LU solver: basic systems")
   - TEST_CASE("LU solver: identity and diagonal")
   - TEST_CASE("LU solver: symmetric positive definite")
   - TEST_CASE("LU solver: stiff systems")
   - TEST_CASE("LU solver: numerically challenging")
   - TEST_CASE("LU solver: known solutions")

4. Build and run the tests:
   cd build
   cmake ..
   make hd_solver_test
   ./hd_solver_test

Notes:
- All test cases use NumPy's solver as a reference
- Epsilon values are adjusted based on the expected numerical difficulty
- Stiff systems include condition numbers ranging from 1e5 to 1e7
- Matrix diagnostic information is provided in comments
""")


if __name__ == "__main__":
    print("/"*80)
    print("// Generated Test Cases for hd_solver.hpp")
    print("// Generated by generate_test_solver.py")
    print("// Reference implementation: NumPy (numpy.linalg.solve)")
    print("/"*80)

    generate_basic_tests()
    generate_identity_tests()
    generate_symmetric_positive_definite_tests()
    generate_known_solution_tests()
    generate_stiff_system_tests()
    generate_tricky_systems()

    print_usage_instructions()
