# Change log

2024/02: add routines for geometric algebra (GA)
2024/04: complete 2d and 3d geometric products for relevant combinations
2024/05: add scripting capabilities for 2d and 3d in lua
2024/07: GA moved to separate project
2024/08: tests for solver added
2025/10: comprehensive test suite expansion and numerical analysis

- Added determinant calculation (hd_determinant.hpp) using LU decomposition
- Comprehensive solver test suite (89 assertions) with stiff systems up to 32×32
- Determinant test suite (24 assertions) with vector and mdspan interfaces
- Python test generators (generate_test_solver.py, generate_test_determinant.py)
- Numerical stability analysis (SOLVER_ANALYSIS.md) confirming near-optimal implementation
- Verified solver handles condition numbers up to 10¹⁰ with appropriate accuracy
- Unified namespace handling (using namespace Kokkos) across solver and determinant
- Developer guide (CLAUDE.md) for future development and AI assistant integration
