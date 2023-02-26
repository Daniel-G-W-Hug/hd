#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

// include functions to be tests
#include "hd_functions.hpp"

TEST_SUITE("fact(n):")
{
    TEST_CASE("fact(n): specific values")
    {
        CHECK(hd::fact(0) == 1.0);
        CHECK(hd::fact(1) == 1.0);
        CHECK(hd::fact(2) == 2.0);
        CHECK(hd::fact(3) == 6.0);
        CHECK(hd::fact(10) == 3628800.);
    }
    TEST_CASE("fact(n): throws on n<0")
    {
        CHECK_THROWS(hd::fact(-1));
    }
}