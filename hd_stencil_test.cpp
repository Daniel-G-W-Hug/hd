#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// include functions to be tested
#include "hd_stencil.hpp"

using namespace hd; // find all functions in hd:: namespace

// Every case pins a stencil with a KNOWN closed form (weights, order, and leading
// truncation coefficient), so a regression in the Taylor-matching system, the
// normalization, or the leading-term detection shows up as a hard numeric mismatch.

TEST_SUITE("stencil_t:")
{
    TEST_CASE("central 3-point first derivative (explicit, order 2)")
    {
        // f'(x0) ~ ( f(x0+h) - f(x0-h) ) / (2h);  leading error  h^2/6 * f'''
        const double x0 = 0.0;
        const double h = 1.0;
        const stencil_t s(x0, stencil_lhs::f1, {x0 - h, x0, x0 + h}, {x0}, {});

        REQUIRE(s.wf0.size() == 3);
        CHECK(s.wf0[0] == doctest::Approx(-1.0 / (2.0 * h)).epsilon(1e-12));
        CHECK(s.wf0[1] == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(s.wf0[2] == doctest::Approx(1.0 / (2.0 * h)).epsilon(1e-12));

        REQUIRE(s.wf1.size() == 1);
        CHECK(s.wf1[0] == doctest::Approx(1.0).epsilon(1e-12)); // normalized lhs

        CHECK(s.order == 2);
        CHECK(s.trunc_err == doctest::Approx(h * h / 6.0).epsilon(1e-12));
    }

    TEST_CASE("central 3-point second derivative (explicit, order 2)")
    {
        // f''(x0) ~ ( f(x0-h) - 2 f(x0) + f(x0+h) ) / h^2;  leading error h^2/12 * f''''
        const double x0 = 0.0;
        const double h = 1.0;
        const stencil_t s(x0, stencil_lhs::f2, {x0 - h, x0, x0 + h}, {}, {x0});

        REQUIRE(s.wf0.size() == 3);
        CHECK(s.wf0[0] == doctest::Approx(1.0 / (h * h)).epsilon(1e-12));
        CHECK(s.wf0[1] == doctest::Approx(-2.0 / (h * h)).epsilon(1e-12));
        CHECK(s.wf0[2] == doctest::Approx(1.0 / (h * h)).epsilon(1e-12));

        REQUIRE(s.wf2.size() == 1);
        CHECK(s.wf2[0] == doctest::Approx(1.0).epsilon(1e-12)); // normalized lhs

        CHECK(s.order == 2);
        CHECK(s.trunc_err == doctest::Approx(h * h / 12.0).epsilon(1e-12));
    }

    TEST_CASE("one-sided 2-point first derivative (explicit, order 1)")
    {
        // f'(x0) ~ ( f(x0+h) - f(x0) ) / h;  leading error  h/2 * f''
        //
        // Pins the leading-term detection: the first non-vanishing residual order is
        // the one that counts (an overwriting scan would report order 2 here).
        const double x0 = 0.0;
        const double h = 1.0;
        const stencil_t s(x0, stencil_lhs::f1, {x0, x0 + h}, {x0}, {});

        REQUIRE(s.wf0.size() == 2);
        CHECK(s.wf0[0] == doctest::Approx(-1.0 / h).epsilon(1e-12));
        CHECK(s.wf0[1] == doctest::Approx(1.0 / h).epsilon(1e-12));

        CHECK(s.order == 1);
        CHECK(s.trunc_err == doctest::Approx(h / 2.0).epsilon(1e-12));
    }

    TEST_CASE("compact 3-point first derivative (implicit Pade scheme, order 4)")
    {
        // (1/6) f'(x0-h) + (2/3) f'(x0) + (1/6) f'(x0+h)
        //     ~ ( f(x0+h) - f(x0-h) ) / (2h)
        //
        // The classic 4th-order compact scheme in the lhs-sum-normalized form.
        // Exercises the lhs (sfact) sign path of the residual computation.
        const double x0 = 0.0;
        const double h = 1.0;
        const stencil_t s(x0, stencil_lhs::f1, {x0 - h, x0 + h}, {x0 - h, x0, x0 + h},
                          {});

        REQUIRE(s.wf0.size() == 2);
        CHECK(s.wf0[0] == doctest::Approx(-1.0 / (2.0 * h)).epsilon(1e-12));
        CHECK(s.wf0[1] == doctest::Approx(1.0 / (2.0 * h)).epsilon(1e-12));

        REQUIRE(s.wf1.size() == 3);
        CHECK(s.wf1[0] == doctest::Approx(1.0 / 6.0).epsilon(1e-12));
        CHECK(s.wf1[1] == doctest::Approx(2.0 / 3.0).epsilon(1e-12));
        CHECK(s.wf1[2] == doctest::Approx(1.0 / 6.0).epsilon(1e-12));
        CHECK(s.wf1[0] + s.wf1[1] + s.wf1[2] == doctest::Approx(1.0).epsilon(1e-12));

        CHECK(s.order == 4);
        CHECK(s.trunc_err == doctest::Approx(-std::pow(h, 4) / 180.0).epsilon(1e-12));
    }

    TEST_CASE("measured convergence rate matches the reported order")
    {
        // apply the one-sided 3-point stencil for f' to f = sin at x0 and measure the
        // convergence rate from two step sizes; it must match the reported order
        const double x0 = 0.3;

        auto fd_error = [&](double h) {
            const stencil_t s(x0, stencil_lhs::f1, {x0, x0 + h, x0 + 2.0 * h}, {x0},
                              {});
            double fd = 0.0;
            for (size_t j = 0; j < s.xf0.size(); ++j) {
                fd += s.wf0[j] * std::sin(s.xf0[j]);
            }
            return std::abs(fd - std::cos(x0));
        };

        // reported order (h large enough for the absolute residual threshold)
        const stencil_t s(x0, stencil_lhs::f1, {x0, x0 + 1.0, x0 + 2.0}, {x0}, {});
        CHECK(s.order == 2);

        const double e1 = fd_error(1.0e-3);
        const double e2 = fd_error(0.5e-3);
        const double rate = std::log2(e1 / e2);
        CHECK(rate == doctest::Approx(2.0).epsilon(0.05));
    }

    TEST_CASE("throws on inconsistent specification")
    {
        // no derivative points at all
        CHECK_THROWS(stencil_t(0.0, stencil_lhs::f1, {0.0, 1.0, 2.0}, {}, {}));
        // lhs f2 requested but no f'' points given
        CHECK_THROWS(stencil_t(0.0, stencil_lhs::f2, {0.0, 1.0}, {0.0}, {}));
        // fewer than 3 points in total
        CHECK_THROWS(stencil_t(0.0, stencil_lhs::f1, {0.0}, {0.0}, {}));
    }
}
