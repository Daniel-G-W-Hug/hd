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

2026/07: stencil fixes and test suite (backport from the ga project)

- Fixed leading-truncation-term detection in hd_stencil.hpp: the residual scan kept
  overwriting order/trunc_err, so the LAST non-vanishing Taylor order was reported
  instead of the first (a one-sided 2-point f' stencil came out as order 2 instead
  of 1); the scan now stops at the first non-vanishing term
- Fixed lhs sign in the truncation residual: lhs terms (compact/implicit schemes) now
  enter with the same sfact sign convention used to build the Taylor-matching matrix
- Initialized order/trunc_err members (were left undefined when every residual term
  stayed below the detection threshold)
- Normalized hd_stencil.hpp includes to sibling style ("hd_functions.hpp" instead of
  "hd/hd_functions.hpp") so the header also compiles inside this repo's own test builds
- New stencil test suite (hd_stencil_test.cpp, 36 assertions): closed-form
  weight/order/truncation gates for central f'/f'', one-sided f' (pins the leading-term
  fix), the classic 4th-order compact Pade f' scheme (pins the sfact fix,
  trunc_err = -h^4/180), a measured-convergence-rate check, and ctor validation
