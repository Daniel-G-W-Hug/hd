# HD Library - Developer Guide for AI Assistants

## Library Overview

**hd** is a header-only C++ library (namespace `hd`) providing numerical computation
utilities and helper functions. All functionality is header-only for ease of use - just
include the relevant `.hpp` file.

**Language Standard**: C++23
**Build System**: CMake + Ninja (for tests only)
**Testing Framework**: doctest

## Project Structure

``` text
hd/
├── *.hpp                          # Header-only library files
├── *_test.cpp                     # Unit tests (not needed for library usage)
├── generate_test_*.py             # Python scripts for generating test cases
├── CMakeLists.txt                 # Build configuration (tests only)
├── README.md                      # Basic usage info
├── CHANGELOG.md                   # Project history
├── SOLVER_ANALYSIS.md            # Numerical analysis of LU solver
└── build/                         # Build artifacts (ignore)
```

## Core Components

### 1. Numerical Solvers & Linear Algebra

#### `hd_solver.hpp` - LU Decomposition Solver

**Status**: Production-ready, optimal implementation from Numerical Recipes

**Functions**:

- `void lu_decomp(mdspan<double, dextents<size_t, 2>> a, mdspan<int, dextents<size_t, 1>> perm)`
- `void lu_backsubs(mdspan<double const, dextents<size_t, 2>> a, mdspan<int const,
  dextents<size_t, 1>> perm, mdspan<double, dextents<size_t, 1>> b)`

**Key Features**:

- Scaled partial pivoting (optimal for general use)
- In-place decomposition
- Handles condition numbers up to ~10¹⁰
- Excellent numerical stability (tested up to 32×32 matrices)

**Test Coverage**: 89 assertions, 7 test cases including stiff systems
**Test Generator**: `generate_test_solver.py`

#### `hd_determinant.hpp` - Matrix Determinant

**Dependencies**: `hd_solver.hpp`

**Functions**:

- `T det(const std::vector<std::vector<T>>& A)` - Vector of vectors interface
- `T det(mdspan<T, Extents, LayoutPolicy, AccessorPolicy> A)` - mdspan interface

**Implementation**: Uses LU decomposition with sign correction from row swaps

**Test Coverage**: 24 assertions, 8 test cases
**Test Generator**: `generate_test_determinant.py`

#### `hd_stencil.hpp` - Finite Difference Stencil Generation

**Dependencies**: `hd_solver.hpp`, `hd_functions.hpp`

**Purpose**: Generate finite difference coefficients for numerical derivatives
**Note**: Uses the LU solver internally

### 2. Mathematical Functions

#### `hd_functions.hpp` - General Math Utilities

**Functions**:

- `double linear_step(double low_x, double high_x, double x)`
- `double smooth_step(double low_x, double high_x, double x)` - Smoothstep interpolation
- `double smoother_step(double low_x, double high_x, double x)` - Smootherstep
- `double log_gamma(double xx)` - Log gamma function
- `double fact(int n)` - Factorial
- `double log_fact(int n)` - Log factorial
- `double bico(int n, int k)` - Binomial coefficient
- `int oo_magnitude(double x, split_t s)` - Order of magnitude

**Special Functions**:

- `template<typename T> constexpr T kronecker(T i, T j)` - Kronecker delta
- `template<typename... T> constexpr int eps(T... indices)` - Levi-Civita permutation
  symbol

**Test Coverage**: Comprehensive tests in `hd_functions_test.cpp`

### 3. Utility Components

#### `hd_error.hpp` - Error Handling

**Features**:

- Source location tracking (`std::source_location`)
- `file_loc_msg()` - Location-aware error messages
- `Error_action` enum for error handling strategies

#### `hd_stop_watch.hpp` - Timing Utilities

**Purpose**: Performance measurement

#### `hd_keypress.hpp` - Console Input

**Purpose**: Non-blocking keyboard input handling

#### `hd_string_trim.hpp` - String Utilities

**Purpose**: String trimming operations

#### `hd_thrdsf_queue.hpp` - Thread-Safe Queue

**Purpose**: Lock-free/thread-safe data structures

#### `hd_thrdsf_stack.hpp` - Thread-Safe Stack

**Purpose**: Lock-free/thread-safe data structures

## Dependencies

### External Libraries

- **mdspan**: C++23 mdspan implementation (Kokkos::mdspan)
  - Location: `../../include/mdspan/include`
  - Used for: Multi-dimensional array views
  - **Important**: Uses `MDSPAN_USE_BRACKET_OPERATOR=1` for `a[i,j]` syntax

- **doctest**: Testing framework (tests only)
  - Install: `brew install doctest`

- **fmt**: Formatting library
  - Install: `brew install fmt`
  - Used in: solver/determinant tests for output

- **date**: Howard Hinnant's date library (commented out)
  - Install: `brew install howard-hinnant-date`

## Coding Conventions

### Namespace Usage

```cpp
// At file scope (before namespace hd):
using namespace Kokkos;  // Makes mdspan less verbose

namespace hd {
    // Library code
}
```

**Rationale**: Consistent with `hd_solver.hpp`, improves readability

### mdspan Syntax

```cpp
// Use bracket operator (enabled via compile flag)
double value = matrix[i, j];  // NOT matrix(i, j)
```

**Compile Flags**:

- `MDSPAN_USE_BRACKET_OPERATOR=1`
- `MDSPAN_CXX_STANDARD=23`

### Error Handling

- Solver: Uses `Solver_error` exception
- Determinant: Catches `Solver_error`, returns 0 for singular matrices
- General: Uses `std::invalid_argument` for precondition violations

## Building and Testing

### Build Commands

```bash
cd build
cmake ..              # Or: cmake .. -GNinja
ninja                 # Or: make
```

### Running Tests

```bash
# All tests
./hd_functions_test
./hd_solver_test
./hd_determinant_test

# Specific test cases
./hd_solver_test -tc="*stiff*"      # Only stiff systems
./hd_solver_test -tc="*32x32*"     # Only 32x32 test
```

### Test Statistics

- **hd_functions_test**: Multiple test suites for math functions
- **hd_solver_test**: 89 assertions, 7 test cases
- **hd_determinant_test**: 24 assertions, 8 test cases

## Test Generation Workflow

### For Solver Tests

1. Edit `generate_test_solver.py` to add new test cases
2. Run: `python3 generate_test_solver.py > output.txt`
3. Copy relevant test cases to `hd_solver_test.cpp`
4. Rebuild and verify tests pass

### For Determinant Tests

1. Edit `generate_test_determinant.py`
2. Run and copy to `hd_determinant_test.cpp`

**Important**: Always verify solutions with NumPy before adding to test suite

## Known Issues and Considerations

### Numerical Accuracy

- LU solver maintains 7-8 digits accuracy for κ ~ 10⁸
- For κ > 10¹⁰, consider iterative refinement (not currently implemented)
- See `SOLVER_ANALYSIS.md` for detailed analysis

### mdspan Version

- Uses Kokkos single-header mdspan (branch: single-header)
- Ensure `MDSPAN_USE_BRACKET_OPERATOR=1` is set in compile flags
- Location: `../../include/mdspan/include`

### Platform-Specific

- Built/tested on macOS with AppleClang
- CMake minimum version: 3.29
- Warnings enabled but not treated as errors (no `-Werror`)

## Future Development Areas

### Potential Enhancements

1. **Iterative Refinement** (for LU solver)
   - Would improve accuracy for κ > 10⁷
   - Cost: ~2x solve time
   - Implementation sketch provided in `SOLVER_ANALYSIS.md`

2. **Cholesky Decomposition**
   - For symmetric positive definite matrices
   - Would be faster than LU for SPD cases
   - Same API pattern as LU solver

3. **QR Decomposition**
   - Better for least squares problems
   - More stable for rank-deficient matrices

4. **Matrix Inversion**
   - Using existing LU decomposition
   - API: `matrix_inverse(A)` returns A⁻¹

5. **Additional Test Coverage**
   - Non-square systems (should throw)
   - Extremely large systems (n > 100)
   - Complex number support

### Documentation Improvements

- Add Doxygen comments to all headers
- Create examples directory with usage patterns
- Document performance characteristics

## Working with This Library

### When Adding New Features

1. Create header file `hd_newfeature.hpp`
2. Add to `namespace hd`
3. Add `using namespace Kokkos;` if using mdspan
4. Create test file `hd_newfeature_test.cpp`
5. Update `CMakeLists.txt` with new test executable
6. Update this file (CLAUDE.md) with new component info

### When Modifying Existing Code

1. Always run full test suite after changes
2. Check numerical accuracy hasn't degraded
3. Verify consistency with coding conventions
4. Update documentation if API changes

### Git Workflow Notes

Current uncommitted changes:

- Modified: `CMakeLists.txt`
- Modified: `hd_functions.hpp`
- Modified: `hd_functions_test.cpp`
- Modified: `hd_determinant.hpp` (namespace cleanup)
- New: `generate_test_determinant.py`
- New: `hd_determinant.hpp`, `hd_determinant_test.cpp`
- New: `generate_test_solver.py`, `hd_solver_test.cpp`
- New: `SOLVER_ANALYSIS.md`, `CLAUDE.md`

## Quick Reference

### Common Tasks

**Add a new test to solver**:

```cpp
SUBCASE("My new test") {
    std::array<double, 9> m_s{...};  // Matrix data
    std::array<double, 3> rhs_s{...}; // RHS
    std::array x_expected{...};       // Solution
    std::array<int, 3> m_perm_s;

    auto m = mdspan<double, extents<size_t, 3, 3>>(m_s.data());
    auto rhs = mdspan<double, extents<size_t, 3>>(rhs_s.data());
    auto m_perm = mdspan<int, extents<size_t, 3>>(m_perm_s.data());

    hd::lu_decomp(m, m_perm);
    hd::lu_backsubs(m, m_perm, rhs);

    for (size_t i = 0; i < 3; ++i) {
        CHECK(std::abs(rhs[i] - x_expected[i]) < 1e-12);
    }
}
```

**Verify with NumPy**:

```python
import numpy as np
A = np.array([[...]])
b = np.array([...])
x = np.linalg.solve(A, b)
print(f"Condition number: {np.linalg.cond(A)}")
print(f"Solution: {x}")
```

## Version History

- **2024/02**: Add geometric algebra routines
- **2024/04**: Complete 2D/3D geometric products
- **2024/05**: Add Lua scripting for 2D/3D
- **2024/07**: GA moved to separate project
- **2024/08**: Tests for solver added
- **2024/10**: Comprehensive determinant & solver tests, numerical analysis, namespace cleanup

---

*This file is intended to help AI assistants quickly understand and work with the hd library. Update as the library evolves.*
