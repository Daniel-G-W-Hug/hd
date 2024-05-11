// author: Daniel Hug, 2024

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <chrono>
#include <tuple>
#include <vector>

// include functions to be tested
#include "hd_ga.hpp"

#include "fmt/chrono.h" // chrono support
#include "fmt/format.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

TEST_SUITE("Geometric Algebra")
{
    using namespace hd::ga;

    TEST_CASE("algebra<2> - 2d_ega")
    {
        fmt::println("");
        fmt::println("algebra<2> - 2d_ega:");
        // 2d euklidean geometric algebra
        const algebra<2> alg;
        CHECK(alg.p() == 2);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 2);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 4);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 3); // == dim_space + 1
        fmt::println("   2d_ega: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   2d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("algebra<3> - 3d_ega")
    {
        fmt::println("");
        fmt::println("algebra<3> - 3d_ega:");
        // 3d euklidean geometric algebra
        const algebra<3> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 3);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 8);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("   3d_ega: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   3d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("algebra<4> 4d_ega")
    {
        fmt::println("");
        fmt::println("algebra<4> - 4d_ega:");
        // 4d euklidean geometric algebra
        const algebra<4> alg;
        CHECK(alg.p() == 4);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   4d_ega: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   4d_ega: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec2d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2d: default init")
    {
        fmt::println("Vec2d: default init");
        Vec2d v;
        // fmt::println("   v = {}", v);
        CHECK(v.x == 0.0);
        CHECK(v.y == 0.0);
    }
    TEST_CASE("Vec2d: with curly braced intializer")
    {
        fmt::println("Vec2d: with curly braced intializer");
        Vec2d v{0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(v.x == 0.0);
        CHECK(v.y == 0.0);
    }
    TEST_CASE("Vec2d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec2d: cp ctor & cp assign incl. type deduction");
        Vec2d v1{1.0, 2.0}; // init with double (type deduction)
        Vec2d v2{v1};       // cp ctor
        Vec2d v3 = v2;      // cp assign
        Vec2d v4 = -v2;     // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(v1.x == 1.0);
        CHECK(v1.y == 2.0);
        CHECK(v2.x == 1.0);
        CHECK(v2.y == 2.0);
        CHECK(v3.x == 1.0);
        CHECK(v3.y == 2.0);
    }


    TEST_CASE("Vec2d: fmt & cout printing")
    {
        fmt::println("Vec2d: fmt & cout printing");

        Vec2d pf{1.0f, 2.0001f};
        Vec2d pd{1.0, 2.0001};

        std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.8f}", pf);

        std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.8f}", pd);

        std::vector<Vec2d<double>> vp1{{1.0, 1.0}, {1.5, 2.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:.e}", fmt::join(vp1, ", "));
    }

    TEST_CASE("Vec2d: comparison float")
    {
        fmt::println("");
        fmt::println("Vec2d: comparison float");

        Vec2d<float> v1f{1.0f, 2.0f};
        Vec2d<float> v2f{2.0f, 4.0f};
        Vec2d<float> v3f{1.0f, 2.0000001f};
        Vec2d<float> v4f{v1f};

        // fmt::println("   v1f = {}", v1f);
        // fmt::println("   v2f = {}", v2f);
        // fmt::println("   v3f = {}", v3f);
        // fmt::println("   v4f = {}", v4f);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<float>::epsilon());

        CHECK(v1f == v4f);           // comparison (equality)
        CHECK(v1f != v2f);           // comparison (inequality)
        CHECK(nrm(v1f) < nrm(v2f));  // comparison (less than)
        CHECK(nrm(v2f) >= nrm(v1f)); // comparison (greater than or equal)
        CHECK(v3f == v1f);           // comparison (eqality)
    }

    TEST_CASE("Vec2d: comparison double")
    {
        fmt::println("Vec2d: comparison double");

        Vec2d<double> v1d{1.0, 2.0};
        Vec2d<double> v2d{2.0, 4.0};
        Vec2d<double> v3d{1.0, 2.0000000000000001};
        Vec2d<double> v4d{v1d};

        // fmt::println("   v1d = {}", v1d);
        // fmt::println("   v2d = {}", v2d);
        // fmt::println("   v3d = {}", v3d);
        // fmt::println("   v4d = {}", v4d);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<double>::epsilon());

        CHECK(v1d == v4d);           // comparison (equality)
        CHECK(v1d != v2d);           // comparison (inequality)
        CHECK(nrm(v1d) < nrm(v2d));  // comparison norm
        CHECK(nrm(v2d) >= nrm(v1d)); // comparison norm
        CHECK(v3d == v1d);           // comparison (eqality)
    }

    TEST_CASE("Vec2d: vector space and linearity tests")
    {
        fmt::println("Vec2d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        Vec2d p0;
        Vec2d p1{1.0, 2.0};
        Vec2d p2{2.0, 4.0};
        Vec2d p3{3.0, 6.0};
        Vec2d p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t. multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t. addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
    }

    TEST_CASE("Vec2d: inner product properties")
    {
        fmt::println("Vec2d: inner product properties");

        double a = 2.35;
        Vec2d u{1.0, 2.0};
        Vec2d v{-0.5, 3.0};
        Vec2d w{3.0, 6.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec2d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2d: operations - norm, inverse, dot")
    {
        fmt::println("Vec2d: operations - norm, inverse, dot");

        Vec2d v1{2.0, 1.0};
        Vec2d v2{unitized(v1)};

        vec2d v3{2.0, 6.0};
        vec2d v4{inv(v3)};

        // fmt::println("v1 = {: .8f}, nrm(v1) = {: .8f}", v1, nrm(v1));
        // fmt::println("v2 = unitized(v1) = {: .8f}, nrm(v2) = {: .8f}", v2, nrm(v2));

        CHECK(std::abs(sq_nrm(v1) - 5.0) < eps);
        CHECK(std::abs(sq_nrm(v2) - 1.0) < eps);
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);
    }

    TEST_CASE("Vec2d: operations - angle")
    {
        fmt::println("Vec2d: operations - angle");

        std::vector<std::tuple<double, Vec2d<double>>> v1;
        std::vector<std::tuple<double, Vec2d<double>>> v2;
        std::vector<std::tuple<double, Vec2d<double>>> v3;

        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2d<double>(std::cos(phi), std::sin(phi));
            v1.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e1_2d, c));
        }
        // fmt::println("");

        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2d<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2));
            v2.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e2_2d, c));
        }
        fmt::println("");

        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2d<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4));
            v3.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e1_2d + e2_2d, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v1) {
            CHECK(std::abs(phi - angle(e1_2d, c)) < eps);
        }
        for (auto const& [phi, c] : v2) {
            CHECK(std::abs(phi - angle(e2_2d, c)) < eps);
        }
        auto ref_vec = unitized(e1_2d + e2_2d);
        for (auto const& [phi, c] : v3) {
            CHECK(std::abs(phi - angle(ref_vec, c)) < eps);
        }

        // auto v = Vec2d<double>{1.0, 0.0};
        // // auto v = Vec2d<double>{1.0, 1.0};
        // for (auto const& [phi, c] : v1) {
        //     auto u1 = v * c;
        //     auto u2 = c * v;
        //     fmt::println("   phi={: .4f}, phi={:> 4.0f}°, c={: .3f},"
        //                  "  u1={: .3f}, u2={: .3f}",
        //                  phi, phi * 180 / pi, c, u1, u2);
        // }
        // fmt::println("");
    }

    TEST_CASE("Vec2d: operations - wedge")
    {
        fmt::println("Vec2d: operations - wedge");

        std::vector<std::tuple<double, Vec2d<double>>> v;

        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2d<double>(std::cos(phi), std::sin(phi));
            v.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e1_2d, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v) {
            CHECK(std::abs(wdg(e1_2d, c) - sin(angle(e1_2d, c))) < eps);
        }
    }

    TEST_CASE("Vec2d: operations - project / reject")
    {
        fmt::println("Vec2d: operations - project / reject");

        vec2d v1{5.0, 1.0};
        vec2d v2{2.0, 2.0};
        vec2d v2u = unitized(v2);

        vec2d v3{project_onto(v1, v2)};
        vec2d v4{reject_from(v1, v2)};
        vec2d v5{v3 + v4};

        vec2d v6{project_onto_unitized(v1, v2u)};
        vec2d v7{reject_from_unitized(v1, v2u)};
        vec2d v8{v6 + v7};

        // fmt::println("v1  = {: .8f}, nrm(v1) = {: .8f}", v1, nrm(v1));
        // fmt::println("v2  = {: .8f}, nrm(v2) = {: .8f}", v2, nrm(v2));
        // fmt::println("v2u = {: .8f}, nrm(v2) = {: .8f}", v2u, nrm(v2u));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .8f}, nrm(v3) = {: .8f}", v3,
        //              nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .8f}, nrm(v4) = {: .8f}", v4,
        //              nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .8f}, nrm(v5) = {: .8f}", v5,
        //              nrm(v5));
        // fmt::println("v6 = project_onto_unitized(v1, v2u) = {: .8f}, nrm(v6) = {:
        // .8f}",
        //              v6, nrm(v6));
        // fmt::println("v7 = reject_from_unitized(v1, v2u)  = {: .8f}, nrm(v7) = {:
        // .8f}",
        //              v7, nrm(v7));
        // fmt::println("v8 = v6 + v7                        = {: .8f}, nrm(v8) = {:
        // .8f}",
        //              v8, nrm(v8));

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);
        CHECK(v6 + v7 == v8);
        CHECK(v8 == v1);

        // checking time required
        //
        // auto start = std::chrono::system_clock::now();
        // for (size_t i = 0; i < 10000000; ++i) {
        //     vec2d v = reject_from(v1, v2);
        // }
        // auto end = std::chrono::system_clock::now();
        // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        // start); fmt::println("The measurement took {}", elapsed);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec2d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2d: default init")
    {
        fmt::println("MVec2d: default init");
        // default initialization
        MVec2d v;
        // fmt::println("   v = {}", v);
        CHECK(v.c0 == 0.0);
        CHECK(v.c1 == 0.0);
        CHECK(v.c2 == 0.0);
        CHECK(v.c3 == 0.0);
    }
    TEST_CASE("MVec2d: with curly braced intializer")
    {
        fmt::println("MVec2d: with curly braced intializer");
        // default initialization
        MVec2d v{0.0, 1.0, 2.0, 3.0};
        // fmt::println("   v = {}", v);
        CHECK(v.c0 == 0.0);
        CHECK(v.c1 == 1.0);
        CHECK(v.c2 == 2.0);
        CHECK(v.c3 == 3.0);
    }

    TEST_CASE("MVec2d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("MVec2d: cp ctor & cp assign incl. type deduction");
        // default initialization
        MVec2d v1{0.0, 1.0, 2.0, 3.0}; // init with double (type deduction)
        MVec2d v2{v1};                 // cp ctor
        MVec2d v3 = v2;                // cp assign
        MVec2d v4 = -v3;               // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(v2.c0 == 0.0);
        CHECK(v2.c1 == 1.0);
        CHECK(v2.c2 == 2.0);
        CHECK(v2.c3 == 3.0);
        CHECK(v3.c0 == 0.0);
        CHECK(v3.c1 == 1.0);
        CHECK(v3.c2 == 2.0);
        CHECK(v3.c3 == 3.0);
    }

    TEST_CASE("MVec2d: fmt & cout printing")
    {
        fmt::println("MVec2d: fmt & cout printing");

        MVec2d pf{1.0f, 2.0001f, 0.0f, 3.0f};
        MVec2d pd{1.0, 2.0001, 0.0, 3.0};

        // std::cout << "   cout: pf = " << pf << std::endl;
        // fmt::println("    fmt: pf = {}", pf);
        // fmt::println("    fmt: pf = {:.8f}", pf);

        // std::cout << "   cout: pd = " << pd << std::endl;
        // fmt::println("    fmt: pd = {}", pd);
        // fmt::println("    fmt: pd = {:.8f}", pd);

        // std::vector<MVec2d<double>> vp1{{1.0, 1.0, 1.0, 2.0}, {0.5, 1.5, 2.0, 2.5}};
        // fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        // fmt::println("    fmt: vp1 = {:.e}", fmt::join(vp1, ", "));

        CHECK(pf == pd);
    }

    TEST_CASE("MVec2d: vector space and linearity tests")
    {
        fmt::println("MVec2d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        MVec2d p0;
        MVec2d p1{0.0, 1.0, 2.0, 3.0};
        MVec2d p2{0.0, 2.0, 4.0, 6.0};
        MVec2d p3{0.0, 3.0, 6.0, 9.0};
        MVec2d p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t. multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t. addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec2d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2d: geometric product tests")
    {
        fmt::println("MVec2d: geometric product tests");

        Vec2d v1{1.0, 2.0};
        Vec2d v2{0.5, 3.0};
        auto d12 = dot(v1, v2);
        auto w12 = wdg(v1, v2);

        MVec2d<double> mv0;
        MVec2d mv1{0.0, 1.0, 2.0, 0.0};
        MVec2d mv2{0.0, 0.5, 3.0, 0.0};
        auto wdp_mv12 = 0.5 * (gpr(mv1, mv2) + gpr(mv2, mv1));
        auto wdm_mv12 = 0.5 * (gpr(mv1, mv2) - gpr(mv2, mv1));

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   dot(v1,v2) = {}", d12);
        // fmt::println("   wdg(v1,v2) = {}", w12);
        // fmt::println("");
        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("   wdp_mv12 = 0.5*(mv1 mv2 + mv2 mv1) = {}", wdp_mv12);
        // fmt::println("   wdm_mv12 = 0.5*(mv1 mv2 - mv2 mv1) = {}", wdm_mv12);
        // fmt::println("");
        // fmt::println("   gr0(wdp_mv12) = {}", gr0(wdp_mv12));
        // fmt::println("   gr1(wdp_mv12) = {}", gr1(wdp_mv12));
        // fmt::println("   gr2(wdp_mv12) = {}", gr2(wdp_mv12));
        // fmt::println("");
        // fmt::println("   gr0(wdm_mv12) = {}", gr0(wdm_mv12));
        // fmt::println("   gr1(wdm_mv12) = {}", gr1(wdm_mv12));
        // fmt::println("   gr2(wdm_mv12) = {}", gr2(wdm_mv12));

        CHECK(d12 == gr0(wdp_mv12));
        CHECK(w12 == gr2(wdm_mv12));
    }

    TEST_CASE(
        "MVec2d geometric product tests - recovering vectors from the geometric product")
    {
        fmt::println("MVec2d: geometric product tests - recovering vectors from the "
                     "geometric product");

        // Two multivectors mv1 and mv2 formed from vectors v1 and v2.
        // (gr0(mv1)==0 && gr1(mv1) != 0 && gr2(mv1)==0 &&
        //  gr0(mv2)==0 && gr1(mv2) != 0 && gr2(mv2)==0 )
        //
        // They are multiplied by the geometric product to form a multivector C
        // C = mv1(v1) mv2(v2) = gpr(mv1,mv2)
        //
        // C contains a scalar part and a bivector part exclusively,
        // the remaining components are zero.
        // (gr0(C) != 0 && gr1(C)==0 && gr2(C) !=0)
        //
        // The scalar part of C represents the parts of v1 and v2
        // that are parallel to each other.
        // The bivector part of C represents the parts of v1 and v2
        // that are perpendicular to each other.
        //
        // multiply C from the right with inv(v2) recovers v1
        // multiply C from the left the the inv(v1) recovers v2

        Vec2d v1{1.0, 2.0};
        Vec2d v2{0.5, 3.0};
        MVec2d mv1{v1};
        MVec2d mv2{v2};

        auto d12 = dot(v1, v2);
        auto w12 = wdg(v1, v2);
        auto mv12 = gpr(mv1, mv2);
        MVec2d mv12a{Scalar2d<double>(d12), PScalar2d<double>(w12)};

        auto v2i = inv(v2);
        auto nv2 = nrm(v2);
        auto nv2i = nrm(v2i);
        auto mv2i{MVec2d(v2i)};
        auto mv1i{MVec2d(inv(v1))};

        // fmt::println("   v1  = {}", v1);
        // fmt::println("   v2  = {}", v2);
        // fmt::println("");
        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("");
        // fmt::println("   dot(v1,v2)   = {}", d12);
        // fmt::println("   wdg(v1,v2)   = {}", w12);
        // fmt::println("   gpr(mv1,mv2) = {}", mv12);
        // fmt::println("   mv12a        = {}", mv12a);
        // fmt::println("");
        // fmt::println("   nv2  = {}", nv2);
        // fmt::println("   v2i  = {}", v2i);
        // fmt::println("   nv2i = {}", nv2i);
        // fmt::println("   mv2i = {}", mv2i);
        // fmt::println("");
        // fmt::println("   gpr(mv12,mv2i) = {}", gpr(mv12, mv2i));
        // fmt::println("   gpr(mv1i,mv12) = {}", gpr(mv1i, mv12));
        // fmt::println("");

        // shorter version w/o intermediate results
        Vec2d a{1.0, 2.0};
        Vec2d b{0.5, 3.0};
        MVec2d C{Scalar2d<double>(dot(a, b)), PScalar2d<double>(wdg(a, b))};
        MVec2d gpr_right = gpr(C, MVec2d{inv(b)});
        MVec2d gpr_left = gpr(MVec2d{inv(a)}, C);

        // fmt::println("   a  = {}", a);
        // fmt::println("   b  = {}", b);
        // fmt::println("   C = gpr(a,b) = {}", C);
        // fmt::println("");
        // fmt::println("   gpr(C,bi) = gpr_right = {}", gpr_right);
        // fmt::println("   gpr(ai,C) = gpr_left = {}", gpr_left);
        // fmt::println("   gr1(gpr_right) = a = {}", gr1(gpr_right));
        // fmt::println("   gr1(gpr_left) = b = {}", gr1(gpr_left));

        CHECK(a == gr1(gpr_right));
        CHECK(b == gr1(gpr_left));
    }

    TEST_CASE("MVec2d: geometric product tests - equivalence tests")
    {
        fmt::println("MVec2d: geometric product tests - equivalence tests");

        Vec2d a{1.0, 2.0};
        Vec2d b{0.5, 3.0};
        MVec2d mva{a};
        MVec2d mvb{b};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        MVec2d ab = gpr(a, b);
        MVec2d abm = gpr(mva, mvb);
        MVec2d abd{Scalar2d<double>(dot_ab), wdg_ab};

        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   ab  = gpr(a,b)                   = {}", ab);
        // fmt::println("   abm = gpr(mva,mvb)               = {}", abm);
        // fmt::println("   abd = MVec2d(dot(a,b), wdg(a,b)) = {}", abd);

        CHECK(ab == abm);
        CHECK(ab == abd);
    }

    TEST_CASE("MVec2d: assignment tests")
    {
        fmt::println("MVec2d: assignment tests");

        Vec2d v1{1.0, 2.0};
        Vec2d v2{0.5, 3.0};
        auto d12 = dot(v1, v2);
        auto w12 = wdg(v1, v2);

        MVec2d<double> mv0;
        MVec2d mv1{0.0, 1.0, 2.0, 0.0};
        MVec2d mv2{0.0, 0.5, 3.0, 0.0};
        MVec2d mv3{v1};
        MVec2d mv4 = v2;

        MVec2d mv5(Scalar2d<double>(5.0));
        MVec2d mv6(PScalar2d<double>(6.0));

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("");
        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("   mv3 = {}", mv3);
        // fmt::println("   mv4 = {}", mv4);
        // fmt::println("   mv5 = {}", mv5);
        // fmt::println("   mv6 = {}", mv6);
        // fmt::println("");
        // fmt::println("   gr1(mv1) = {}", gr1(mv1));
        // fmt::println("   gr1(mv2) = {}", gr1(mv2));
        // fmt::println("   gr1(mv3) = {}", gr1(mv3));
        // fmt::println("   gr1(mv3) = {}", gr1(mv4));

        CHECK(gr1(mv1) == v1);
        CHECK(gr1(mv2) == v2);
        CHECK(gr1(mv3) == v1);
        CHECK(gr1(mv4) == v2);
        CHECK(mv1 == mv3);
        CHECK(mv4 == mv2);
    }

    TEST_CASE("MVec2d: modelling complex numbers")
    {
        fmt::println("MVec2d: modelling complex numbers");


        Vec2d v1{1.0, -1.0};
        MVec2d v1m{v1}; // full 2d multivector

        // multiplying with e1 from the left should make it a complex number
        // i.e. a multivector with a scalar (=Re) and a bivector part (=Im)
        // (for test purposes here, the even subalgebra would be sufficient)
        auto vc = gpr(e1_2d, v1);
        auto vcm = gpr(e1m_2d, v1m); // full gpr

        // multiplying with I2 from the right should rotate by +90°
        auto vr = gpr(vc, I_2d);
        auto vrm = gpr(vcm, Im_2d); // full gpr

        // multiplying with I2 from the left should rotate by -90°
        auto vl = gpr(I_2d, vc);
        auto vlm = gpr(Im_2d, vcm); // full gpr

        // defining a complex number in all three forms
        auto u = Vec2d{1.0, 0.0};
        auto v = Vec2d{std::cos(pi / 6.0), std::sin(pi / 6.0)}; // unit vec +30%
        auto angle_uv = angle(u, v);

        auto uv = gpr(u, v); // complex number with real part and bivector part
        auto a = gr0(uv);
        auto b = gr2(uv);
        auto r = sqrt(a * a + b * b);

        // fmt::println("   I_2d             = {}", I_2d);
        // fmt::println("   Im_2d            = {}", Im_2d);
        // fmt::println("   grp(I_2d,I_2d)   = {}", gpr(I_2d, I_2d));
        // fmt::println("   grp(Im_2d,Im_2d) = {}", gpr(Im_2d, Im_2d));
        // fmt::println("");
        // fmt::println("   e1_2d  = {}", e1_2d);
        // fmt::println("   e1m_2d = {}", e1m_2d);
        // fmt::println("   e2_2d  = {}", e2_2d);
        // fmt::println("   e2m_2d = {}", e2m_2d);
        // fmt::println("");
        // fmt::println("   vc   = {}", vc);
        // fmt::println("   vcm  = {}", vcm);
        // fmt::println("   vr   = {}", vr);
        // fmt::println("   vrm  = {}", vrm);
        // fmt::println("   vl   = {}", vl);
        // fmt::println("   vlm  = {}", vlm);
        // fmt::println("");
        // fmt::println("   v1            = {}", v1);
        // fmt::println("   gpr(v1,I_2d)  = {}", gpr(v1, I_2d));
        // fmt::println("   gpr(I_2d,v1)  = {}", gpr(I_2d, v1));
        // fmt::println("");
        // fmt::println("   u        = {}", u);
        // fmt::println("   v        = {}", v);
        // fmt::println("   angle_uv = {:.3}°", angle_uv * 180 / pi);
        // fmt::println("");
        // fmt::println("   uv                  = {}", uv);
        // fmt::println("   a = gr0(uv)         = {}", a);
        // fmt::println("   b = gr2(uv)         = {}", b);
        // fmt::println("   r = sqrt(a^2 + b^2) = {}", r);
        // fmt::println("   r exp(angle_uv) = {}", r * exp(PScalar2d<double>(angle_uv)));
        // HINT:declaring angle_uv a PScalar2d makes it a bivector angle,
        // i.e. a multiple of the bivector I_2d
        // ATTENTION: if you don't declare it as such, the normal exponential function
        //            will be called, resulting in a scalar result!

        CHECK(gr0(vc) == gr0(vcm));
        CHECK(gr2(vc) == gr2(vcm));
        CHECK(gr0(vr) == gr0(vrm));
        CHECK(gr2(vr) == gr2(vrm));
        CHECK(gr0(vl) == gr0(vlm));
        CHECK(gr2(vl) == gr2(vlm));
        CHECK(v1.x == gpr(v1, I_2d).y); // rotation +90°
        CHECK(v1.y == -gpr(v1, I_2d).x);
        CHECK(v1.x == -gpr(I_2d, v1).y); // rotation -90°
        CHECK(v1.y == gpr(I_2d, v1).x);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec2d_E<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2d_E: modelling complex numbers - basics")
    {
        fmt::println("MVec2d_E: modelling complex numbers - basics");

        // defining a complex number in all three forms as multivector
        auto u = Vec2d{1.0, 0.0};
        auto v = Vec2d{std::cos(pi / 6.0), std::sin(pi / 6.0)}; // unit vec +30°

        auto angle_uv = angle(u, v);

        auto uv = gpr(u, v); // complex number with real part and bivector part
        auto v2 = exp(I_2d, angle_uv);
        auto re = gr0(uv);
        auto im = gr2(uv);
        auto r = sqrt(re * re + im * im);

        MVec2d_E a{1.0, 0.0};
        MVec2d_E b{1.0, 1.0};
        MVec2d_E c{a + b};
        MVec2d_E d{a - b};
        MVec2d_E e{2.0 * b};
        MVec2d_E f{b * 2.0};
        MVec2d_E g{-e};
        MVec2d_E h{0.0, 1.0};
        MVec2d_E as = gpr(a, a);
        MVec2d_E hs = gpr(h, h);

        MVec2d_E i = gpr(b, c);
        MVec2d_E j = b * c; // defined operator*(MVec2d_E,MVec2d_E) as gpr(a,b)
        auto k = I_2d;
        MVec2d_E l = exp(I_2d, pi / 2);
        MVec2d_E m = Im_2d_E;
        MVec2d n = Im_2d;
        // fmt::println("   Multivector form of complex numbers:");
        // fmt::println("   u                     = {}", u);
        // fmt::println("   v                     = {}", v);
        // fmt::println("   angle(u,v)            = {:.3}°", angle_uv * 180 / pi);
        // fmt::println("   uv = gpr(u,v)         = {}", uv);
        // fmt::println("   re = gr0(uv)          = {}", re);
        // fmt::println("   im = gr2(uv)          = {}", im);
        // fmt::println("   r = sqrt(re^2 + im^2) = {}", r);
        // fmt::println("");
        // fmt::println("   Using the even subalgebra only (std form of complex
        // numbers):");
        // // declaring angle_uv a PScalar2d makes it a bivector angle,
        // // i.e. a multiple of the bivector I_2d
        // // ATTENTION: if you don't declare it as such, the normal exponential function
        // //            will be called, resulting in a scalar result!
        // fmt::println("   v2=exp(angle_uv) = {}", v2);
        // fmt::println("");
        // fmt::println("   a         = {}", a);
        // fmt::println("   b         = {}", b);
        // fmt::println("   c = a+b   = {}", c);
        // fmt::println("   d = a-b   = {}", d);
        // fmt::println("   e = 2.0*b = {}", e);
        // fmt::println("   f = b*2.0 = {}", f);
        // fmt::println("   g = -e    = {}", g);
        // fmt::println("");
        // fmt::println("   h =            = {}", h);
        // fmt::println("   as = gpr(a,a)  = {}", as);
        // fmt::println("   hs = gpr(h,h)  = {}", hs);
        // fmt::println("   b*h = gpr(b,h) = {}", gpr(b, h));
        // fmt::println("   h*b = gpr(h,b) = {}", gpr(h, b));
        // fmt::println("   b*h = b*h      = {}", b * h);
        // fmt::println("   h*b = h*b      = {}", h * b);
        // fmt::println("");
        // fmt::println("   i = gpr(b,c) = {}", i);
        // fmt::println("   j = b*c      = {}", j);
        // fmt::println("");
        // fmt::println("   k = I_2d                         = {}", k);
        // fmt::println("   l = exp(PScalar2d<double>(pi/2)) = {:.3f}", l);
        // fmt::println("   m = Im_2d_E                        = {}", m);
        // fmt::println("   n = Im_2d                        = {}", n);

        CHECK(c == a + b);
        CHECK(d == a - b);
        CHECK(e == 2.0 * b);
        CHECK(f == b * 2.0);
        CHECK(g == -e);
        CHECK(as == a);
        CHECK(hs == MVec2d_E(-1.0, 0.0));
        CHECK(v.x == v2.c0);
        CHECK(v.y == v2.c1);
        CHECK(gpr(b, h) ==
              gpr(h, b)); // in 2d the pseudoscalar commutes commutes with complex numbers
        CHECK(i == j);
        CHECK(l == m);
        CHECK(rev(b + c) == rev(b) + rev(c));
        CHECK(rev(b * c) == rev(b) * rev(c));
        CHECK(nrm(b * c) == nrm(b) * nrm(c));
        CHECK(b * c == c * b);

        CHECK(sq_nrm(MVec2d_E{1.0, 1.0}) == 2.0);
        CHECK(nrm(MVec2d_E{1.0, 1.0}) == std::sqrt(2.0));
        CHECK(rev(MVec2d_E{1.0, 1.0}) == MVec2d_E{1.0, -1.0});
        CHECK(std::abs(nrm(unitized(MVec2d_E{1.0, 1.0})) - 1.0) < eps);

        CHECK(MVec2d_E{-1.0, 1.0} * inv(MVec2d_E{-1.0, 1.0}) == MVec2d_E{1.0, 0.0});
        CHECK(gr0(MVec2d_E{-1.0, 1.0} * rev(MVec2d_E{-1.0, 1.0})) ==
              sq_nrm(MVec2d_E{-1.0, 1.0}));
        CHECK(std::abs(gr2(MVec2d_E{-1.0, 1.0} * rev(MVec2d_E{-1.0, 1.0}))) < eps);

        CHECK(angle(MVec2d_E{1.0, 0.0}) == 0);
        CHECK(angle(MVec2d_E{1.0, 1.0}) == pi / 4.0);
        CHECK(angle(MVec2d_E{0.0, 1.0}) == pi / 2.0);
        CHECK(angle(MVec2d_E{-1.0, 1.0}) == pi * 3.0 / 4.0);
        CHECK(angle(MVec2d_E{-1.0, 0.0}) == pi);
        CHECK(angle(MVec2d_E{1.0, -1.0}) == -pi / 4.0);
        CHECK(angle(MVec2d_E{0.0, -1.0}) == -pi / 2.0);
        CHECK(angle(MVec2d_E{-1.0, -1.0}) == -pi * 3.0 / 4.0);

        CHECK(Vec2d{1.0, 0.0} * Vec2d{1.1, 1.1} ==
              rev(Vec2d{1.1, 1.1} * Vec2d{1.0, 0.0}));
        CHECK(exp(I_2d, pi / 4) == rev(exp(I_2d, -pi / 4)));
        CHECK(exp(I_2d, -angle_uv) * u == u * exp(I_2d, angle_uv)); // 2d rotation direct
        CHECK(exp(I_2d, -angle_uv) * u == v);
        CHECK(rotate(u, rotor(I_2d, angle_uv)) == v); // 2d rotation with double product
                                                      // completely as in the 3d case
                                                      // more effort computationally,
                                                      // but independent of dimension of
                                                      // space
    }

    TEST_CASE("MVec2d_E: modelling complex numbers - products")
    {
        fmt::println("MVec2d_E: modelling complex numbers - products");

        // std::vector<std::tuple<double, MVec2d_E<double>>> c_v;
        // for (int i = -12; i <= 12; ++i) {
        //     double phi = i * pi / 12;
        //     MVec2d_E c = exp(PScalar2d<double>(phi));
        //     c_v.push_back(std::make_tuple(phi, c));
        //     fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
        //                  " angle={: .4f}",
        //                  i, phi, phi * 180 / pi, c, angle(c));
        // }
        // fmt::println("");

        // auto v = Vec2d<double>{1.0, 0.0};
        // // auto v = Vec2d<double>{1.0, 1.0};
        // for (auto const& [phi, c] : c_v) {
        //     auto u1 = v * c;
        //     auto u2 = c * v;
        //     fmt::println("   phi={: .4f}, phi={:> 4.0f}°, c={: .3f},"
        //                  "  u1={: .3f}, u2={: .3f}",
        //                  phi, phi * 180 / pi, c, u1, u2);
        // }
        // fmt::println("");


        CHECK(MVec2d_E{2.0, 3.0} * MVec2d{-1.0, 1.5, -2.0, -3.0} ==
              MVec2d{2.0, 0.0, 0.0, 3.0} * MVec2d{-1.0, 1.5, -2.0, -3.0});
        CHECK(MVec2d_E{2.0, 3.0} * Vec2d{1.5, -2.0} ==
              gr1(MVec2d{2.0, 0.0, 0.0, 3.0} * MVec2d{0.0, 1.5, -2.0, 0.0}));

        CHECK(gr0(Vec2d{1.5, -2.0} * Vec2d{2.0, 3.0}) ==
              gr0(MVec2d{0.0, 1.5, -2.0, 0.0} * MVec2d{0.0, 2.0, 3.0, 0.0}));
        CHECK(gr2(Vec2d{1.5, -2.0} * Vec2d{2.0, 3.0}) ==
              gr2(MVec2d{0.0, 1.5, -2.0, 0.0} * MVec2d{0.0, 2.0, 3.0, 0.0}));

        // multiply from left
        CHECK(PScalar2d<double>(1.5) * MVec2d{-1.0, 1.5, -2.0, -3.0} ==
              MVec2d{0.0, 0.0, 0.0, 1.5} * MVec2d{-1.0, 1.5, -2.0, -3.0});

        CHECK(MVec2d(PScalar2d<double>(1.5) * MVec2d_E{-1.0, -3.0}) ==
              MVec2d{0.0, 0.0, 0.0, 1.5} * MVec2d{-1.0, 0.0, 0.0, -3.0});

        CHECK(MVec2d(PScalar2d<double>(1.5) * Vec2d{-1.0, -3.0}) ==
              MVec2d{0.0, 0.0, 0.0, 1.5} * MVec2d{0.0, -1.0, -3.0, 0.0});

        // multiply from right
        CHECK(MVec2d{-1.0, 1.5, -2.0, -3.0} * PScalar2d<double>(1.5) ==
              MVec2d{-1.0, 1.5, -2.0, -3.0} * MVec2d{0.0, 0.0, 0.0, 1.5});

        CHECK(MVec2d_E{-1.0, -3.0} * MVec2d(PScalar2d<double>(1.5)) ==
              MVec2d{-1.0, 0.0, 0.0, -3.0} * MVec2d{0.0, 0.0, 0.0, 1.5});

        CHECK(MVec2d(Vec2d{-1.0, -3.0} * PScalar2d<double>(1.5)) ==
              MVec2d{0.0, -1.0, -3.0, 0.0} * MVec2d{0.0, 0.0, 0.0, 1.5});

        // two bivectors
        CHECK(MVec2d(Scalar2d<double>(PScalar2d<double>(1.5) * PScalar2d<double>(3.0))) ==
              MVec2d{0.0, 0.0, 0.0, 1.5} * MVec2d{0.0, 0.0, 0.0, 3.0});

        // MVec2d_E tests multiply from left
        CHECK(MVec2d_E{-1.0, -3.0} * MVec2d(-1.0, 1.5, -2.0, -3.0) ==
              MVec2d{-1.0, 0.0, 0.0, -3.0} * MVec2d{-1.0, 1.5, -2.0, -3.0});

        CHECK(MVec2d(MVec2d_E{-1.0, -3.0} * Vec2d(1.5, -2.0)) ==
              MVec2d{-1.0, 0.0, 0.0, -3.0} * MVec2d{0.0, 1.5, -2.0, 0.0});

        // MVec2d_E tests multiply from right
        CHECK(MVec2d(-1.0, 1.5, -2.0, -3.0) * MVec2d_E{-1.0, -3.0} ==
              MVec2d{-1.0, 1.5, -2.0, -3.0} * MVec2d{-1.0, 0.0, 0.0, -3.0});

        CHECK(MVec2d(Vec2d(1.5, -2.0) * MVec2d_E{-1.0, -3.0}) ==
              MVec2d{0.0, 1.5, -2.0, 0.0} * MVec2d{-1.0, 0.0, 0.0, -3.0});

        // multiply two MVec2d_E
        CHECK(MVec2d(MVec2d_E(-3.0, 2.0) * MVec2d_E{-1.0, -3.0}) ==
              MVec2d{-3.0, 0.0, 0.0, 2.0} * MVec2d{-1.0, 0.0, 0.0, -3.0});
    }


    TEST_CASE("MVec2d: dualization")
    {
        fmt::println("MVec2d: dualization");


        Vec2d v{1.0, 2.0};                    // 2d vector
        MVec2d vm{10.0, 1.0, 2.0, 30.0};      // full 2d multivector
        MVec2d vm_even{10.0, 0.0, 0.0, 30.0}; // full 2d multivector - even content
        MVec2d_E vm_E{10.0, 30.0};            // even grade 2d multivector

        // dual by left multiplication with Im_2d
        // as defined in Doran/Lasenby "GA for physicists"

        auto vm_dual_manual = Im_2d * vm;
        auto vm_dual = dual(vm);

        auto vm_dual_even_manual = Im_2d * vm_even;
        auto vm_dual_even = dual(vm_even);

        auto vm_dual_manual_E = Im_2d_E * vm_E;
        auto vm_dual_E = dual(vm_E);

        auto v_dual_manual = I_2d * v;
        auto v_dual = dual(v);

        // fmt::println("   I_2d    = {}", I_2d);
        // fmt::println("   Im_2d   = {}", Im_2d);
        // fmt::println("   Im_2d_E = {}", Im_2d_E);
        // fmt::println("");
        // fmt::println("   vm            = {}", vm);
        // fmt::println("   Im_2d*vm      = {}", vm_dual_manual);
        // fmt::println("   dual(vm)      = {}", vm_dual);
        // fmt::println("");
        // fmt::println("   vm_even       = {}", vm_even);
        // fmt::println("   Im_2d*vm_even = {}", vm_dual_even_manual);
        // fmt::println("   dual(vm_even) = {}", vm_dual_even);
        // fmt::println("");
        // fmt::println("   vm_E          = {}", vm_E);
        // fmt::println("   Im_2d_E*vm_E  = {}", vm_dual_manual_E);
        // fmt::println("   dual(vm_E)    = {}", vm_dual_E);
        // fmt::println("");
        // fmt::println("   v             = {}", v);
        // fmt::println("   I_2d * v      = {}", v_dual_manual);
        // fmt::println("   dual(v)       = {}", v_dual);

        CHECK(vm_dual == vm_dual_manual);
        CHECK(vm_dual_even == vm_dual_even_manual);
        CHECK(vm_dual_E == vm_dual_manual_E);
        CHECK(v_dual == v_dual_manual);
        CHECK(dual(Scalar2d<double>(5)) == PScalar2d<double>(5));
        CHECK(dual(PScalar2d<double>(5)) == Scalar2d<double>(-5));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3d: default init")
    {
        fmt::println("Vec3d: default init");
        Vec3d v;
        // fmt::println("   v = {}", v);
        CHECK(v.x == 0.0);
        CHECK(v.y == 0.0);
        CHECK(v.z == 0.0);
    }
    TEST_CASE("Vec3d: with curly braced intializer")
    {
        fmt::println("Vec3d: with curly braced intializer");
        Vec3d v{0.0, 0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(v.x == 0.0);
        CHECK(v.y == 0.0);
        CHECK(v.z == 0.0);
    }
    TEST_CASE("Vec3d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec3d: cp ctor & cp assign incl. type deduction");
        Vec3d v1{1.0, 2.0, 3.0}; // init with double (type deduction)
        Vec3d v2{v1};            // cp ctor
        Vec3d v3 = v2;           // cp assign
        Vec3d v4 = -v2;          // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(v1.x == 1.0);
        CHECK(v1.y == 2.0);
        CHECK(v1.z == 3.0);
        CHECK(v2.x == 1.0);
        CHECK(v2.y == 2.0);
        CHECK(v2.z == 3.0);
        CHECK(v3.x == 1.0);
        CHECK(v3.y == 2.0);
        CHECK(v3.z == 3.0);
    }


    TEST_CASE("Vec3d: fmt & cout printing")
    {
        fmt::println("Vec3d: fmt & cout printing");

        Vec3d pf{1.0f, 2.0001f, 3.0f};
        Vec3d pd{1.0, 2.0001, 3.0};

        std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.8f}", pf);

        std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.8f}", pd);

        std::vector<Vec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:.e}", fmt::join(vp1, ", "));
    }

    TEST_CASE("Vec3d: comparison float")
    {
        fmt::println("");
        fmt::println("Vec3d: comparison float");

        Vec3d<float> v1f{1.0f, 2.0f, 3.0f};
        Vec3d<float> v2f{2.0f, 4.0f, 3.0f};
        Vec3d<float> v3f{1.0f, 2.0000001f, 3.0f};
        Vec3d<float> v4f{v1f};

        // fmt::println("   v1f = {}", v1f);
        // fmt::println("   v2f = {}", v2f);
        // fmt::println("   v3f = {}", v3f);
        // fmt::println("   v4f = {}", v4f);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<float>::epsilon());

        CHECK(v1f == v4f);           // comparison (equality)
        CHECK(v1f != v2f);           // comparison (inequality)
        CHECK(nrm(v1f) < nrm(v2f));  // comparison (less than)
        CHECK(nrm(v2f) >= nrm(v1f)); // comparison (greater than or equal)
        CHECK(v3f == v1f);           // comparison (eqality)
    }

    TEST_CASE("Vec3d: comparison double")
    {
        fmt::println("Vec3d: comparison double");

        Vec3d<double> v1d{1.0, 2.0, 3.0};
        Vec3d<double> v2d{2.0, 4.0, 3.0};
        Vec3d<double> v3d{1.0, 2.0000000000000001, 3.0};
        Vec3d<double> v4d{v1d};

        // fmt::println("   v1d = {}", v1d);
        // fmt::println("   v2d = {}", v2d);
        // fmt::println("   v3d = {}", v3d);
        // fmt::println("   v4d = {}", v4d);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<double>::epsilon());

        CHECK(v1d == v4d);           // comparison (equality)
        CHECK(v1d != v2d);           // comparison (inequality)
        CHECK(nrm(v1d) < nrm(v2d));  // comparison norm
        CHECK(nrm(v2d) >= nrm(v1d)); // comparison norm
        CHECK(v3d == v1d);           // comparison (eqality)
    }

    TEST_CASE("Vec3d: vector space and linearity tests")
    {
        fmt::println("Vec3d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        Vec3d p0;
        Vec3d p1{1.0, 2.0, 3.0};
        Vec3d p2{2.0, 4.0, 6.0};
        Vec3d p3{3.0, 6.0, 9.0};
        Vec3d p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t. multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t.addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
    }

    TEST_CASE("Vec3d: inner product properties")
    {
        fmt::println("Vec3d: inner product properties");

        double a = 2.35;
        Vec3d u{1.0, 2.0, 1.0};
        Vec3d v{-0.5, 3.0, 0.5};
        Vec3d w{3.0, 6.0, -3.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3d: operations - norm, inverse, dot")
    {
        fmt::println("Vec3d: operations - norm, inverse, dot");

        Vec3d<float> v1{2.0, 1.0, 2.0};
        Vec3d v2{unitized(v1)};

        Vec3d v3{2.0, 6.0, -4.0};
        Vec3d v4{inv(v3)};

        // fmt::println("v1 = {: .8f}, nrm(v1) = {: .8f}", v1, nrm(v1));
        // fmt::println("v2 = unitized(v1) = {: .8f}, nrm(v2) = {: .8f}", v2, nrm(v2));
        // fmt::println("v3 = {: .8f}, nrm(v1) = {: .8f}", v3, nrm(v3));
        // fmt::println(
        //     "v4 = inv(v3) = {: .8f}, nrm(v3) = {: .8f}, nrm(v3)*nrm(v4) = {: .8f}", v4,
        //     nrm(v4), nrm(v3) * nrm(v4));

        CHECK(std::abs(sq_nrm(v1) - 9.0) < eps);
        CHECK(std::abs(sq_nrm(v2) - 1.0) < eps);
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);
    }

    TEST_CASE("Vec3d: operations - angle I")
    {
        fmt::println("Vec3d: operations - angle");

        Vec3d v1{1.0, 0.0, 0.0};
        Vec3d v2{unitized(Vec3d(1.0, 1.0, 0.0))};
        Vec3d v3{0.0, 1.0, 0.0};
        Vec3d v4{unitized(Vec3d(-1.0, 1.0, 0.0))};
        Vec3d v5{-1.0, 0.0, 0.0};
        Vec3d v6{unitized(Vec3d(-1.0, -1.0, 0.0))};
        Vec3d v7{0.0, -1.0, 0.0};
        Vec3d v8{unitized(Vec3d(1.0, -1.0, 0.0))};

        // fmt::println("v1 = {: .8f}, nrm(v1) = {:.8f}, angle(v1,v1) = {:.8f}, {:.8f}",
        // v1,
        //              nrm(v1), angle(v1, v1), angle(v1, v1) / pi);
        // fmt::println("v2 = {: .8f}, nrm(v2) = {:.8f}, angle(v1,v2) = {:.8f}, {:.8f}",
        // v2,
        //              nrm(v2), angle(v1, v2), angle(v1, v2) / pi);
        // fmt::println("v3 = {: .8f}, nrm(v3) = {:.8f}, angle(v1,v3) = {:.8f}, {:.8f}",
        // v3,
        //              nrm(v3), angle(v1, v3), angle(v1, v3) / pi);
        // fmt::println("v4 = {: .8f}, nrm(v4) = {:.8f}, angle(v1,v4) = {:.8f}, {:.8f}",
        // v4,
        //              nrm(v4), angle(v1, v4), angle(v1, v4) / pi);
        // fmt::println("v5 = {: .8f}, nrm(v5) = {:.8f}, angle(v1,v5) = {:.8f}, {:.8f}",
        // v5,
        //              nrm(v5), angle(v1, v5), angle(v1, v5) / pi);
        // fmt::println("v6 = {: .8f}, nrm(v6) = {:.8f}, angle(v1,v6) = {:.8f}, {:.8f}",
        // v6,
        //              nrm(v6), angle(v1, v6), angle(v1, v6) / pi);
        // fmt::println("v7 = {: .8f}, nrm(v7) = {:.8f}, angle(v1,v7) = {:.8f}, {:.8f}",
        // v7,
        //              nrm(v7), angle(v1, v7), angle(v1, v7) / pi);
        // fmt::println("v8 = {: .8f}, nrm(v8) = {:.8f}, angle(v1,v8) = {:.8f}, {:.8f}",
        // v8,
        //              nrm(v8), angle(v1, v8), angle(v1, v8) / pi);

        CHECK(std::abs(angle(v1, v1) - 0.0) < eps);
        CHECK(std::abs(angle(v1, v2) - pi * 0.25) < eps);
        CHECK(std::abs(angle(v1, v3) - pi * 0.5) < eps);
        CHECK(std::abs(angle(v1, v4) - pi * 0.75) < eps);
        CHECK(std::abs(angle(v1, v5) - pi) < eps);
    }

    TEST_CASE("Vec3d: operations - angle II")
    {
        fmt::println("Vec3d: operations - angle II");

        std::vector<std::tuple<double, Vec3d<double>>> v1;
        std::vector<std::tuple<double, Vec3d<double>>> v2;
        std::vector<std::tuple<double, Vec3d<double>>> v3;

        // only positive angles are easy to implement vs. the 2d case

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi), std::sin(phi), 0.0);
            v1.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e1_3d, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2), 0.0);
            v2.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e2_3d, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4), 0.0);
            v3.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e1_3d + e2_3d, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v1) {
            CHECK(std::abs(phi - angle(e1_3d, c)) < eps);
        }
        for (auto const& [phi, c] : v2) {
            CHECK(std::abs(phi - angle(e2_3d, c)) < eps);
        }
        auto ref_vec = unitized(e1_3d + e2_3d);
        for (auto const& [phi, c] : v3) {
            CHECK(std::abs(phi - angle(ref_vec, c)) < eps);
        }
    }

    TEST_CASE("Vec3d: operations - wedge")
    {
        fmt::println("Vec3d: operations - wedge");

        Vec3d v1{1.0, 0.0, 0.0};
        Vec3d v2{unitized(Vec3d(1.0, 1.0, 0.0))};
        Vec3d v3{0.0, 1.0, 0.0};
        Vec3d v4{unitized(Vec3d(-1.0, 1.0, 0.0))};
        Vec3d v5{-1.0, 0.0, 0.0};
        Vec3d v6{unitized(Vec3d(-1.0, -1.0, 0.0))};
        Vec3d v7{0.0, -1.0, 0.0};
        Vec3d v8{unitized(Vec3d(1.0, -1.0, 0.0))};

        // fmt::println("v1 = {: .8f}, wdg(v1,v1) = {: .8f}, "
        //              "angle = {: .8f}",
        //              v1, wdg(v1, v1), angle(v1, v1));
        // fmt::println("v2 = {: .8f}, wdg(v1,v2) = {: .8f}, "
        //              "angle = {: .8f}",
        //              v2, wdg(v1, v2), angle(v1, v2));
        // fmt::println("v3 = {: .8f}, wdg(v1,v3) = {: .8f}, "
        //              "angle = {: .8f}",
        //              v3, wdg(v1, v3), angle(v1, v3));
        // fmt::println("v4 = {: .8f}, wdg(v1,v4) = {: .8f}, "
        //              "angle = {: .8f}",
        //              v4, wdg(v1, v4), angle(v1, v4));
        // fmt::println("v5 = {: .8f}, wdg(v1,v5) = {: .8f}, "
        //              "angle = {: .8f}",
        //              v5, wdg(v1, v5), angle(v1, v5));
        // fmt::println("v6 = {: .8f}, wdg(v1,v6) = {: .8f}, "
        //              "angle = {: .8f}",
        //              v6, wdg(v1, v6), angle(v1, v6));
        // fmt::println("v7 = {: .8f}, wdg(v1,v7) = {: .8f}, "
        //              "angle = {: .8f}",
        //              v7, wdg(v1, v7), angle(v1, v7));
        // fmt::println("v8 = {: .8f}, wdg(v1,v8) = {: .8f}, "
        //              "angle = {: .8f}",
        //              v8, wdg(v1, v8), angle(v1, v8));

        CHECK(std::abs(nrm(wdg(v1, v1)) - sin(angle(v1, v1))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v2)) - sin(angle(v1, v2))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v3)) - sin(angle(v1, v3))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v4)) - sin(angle(v1, v4))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v5)) - sin(angle(v1, v5))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v6)) - sin(angle(v1, v6))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v7)) - sin(angle(v1, v7))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v8)) - sin(angle(v1, v8))) < eps);
    }

    TEST_CASE("Vec3d: operations - project / reject (vector - vector)")
    {
        fmt::println("Vec3d: operations - project / reject (vector - vector)");

        Vec3d v1{5.0, 1.0, 1.0};
        Vec3d v2{2.0, 2.0, 1.0};
        Vec3d v2u = unitized(v2);

        Vec3d v3{project_onto(v1, v2)};
        Vec3d v4{reject_from(v1, v2)};
        Vec3d v5{v3 + v4};

        Vec3d v6{project_onto_unitized(v1, v2u)};
        Vec3d v7{reject_from_unitized(v1, v2u)};
        Vec3d v8{v6 + v7};

        // fmt::println("v1  = {: .8f}, nrm(v1) = {: .8f}", v1, nrm(v1));
        // fmt::println("v2  = {: .8f}, nrm(v2) = {: .8f}", v2, nrm(v2));
        // fmt::println("v2u = {: .8f}, nrm(v2) = {: .8f}", v2u, nrm(v2u));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .8f}, nrm(v3) = {: .8f}", v3,
        //              nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .8f}, nrm(v4) = {: .8f}", v4,
        //              nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .8f}, nrm(v5) = {: .8f}", v5,
        //              nrm(v5));
        // fmt::println("v6 = project_onto_unitized(v1, v2u) = {: .8f}, nrm(v6) = {: .8f}
        // ",
        //              v6, nrm(v6));
        // fmt::println("v7 = reject_from_unitized(v1, v2u)  = {: .8f}, nrm(v7) = {: .8f}
        // ",
        //              v7, nrm(v7));
        // fmt::println("v8 = v6 + v7                        = {: .8f}, nrm(v8) = {: .8f}
        // ",
        //              v8, nrm(v8));

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);
        CHECK(v6 + v7 == v8);
        CHECK(v8 == v1);

        // checking time required
        //
        // auto start = std::chrono::system_clock::now();
        // for (size_t i = 0; i < 10000000; ++i) {
        //     Vec3d v = reject_from(v1, v2);
        // }
        // auto end = std::chrono::system_clock::now();
        // auto elapsed =
        // std::chrono::duration_cast<std::chrono::milliseconds>(end
        // -
        // start); fmt::println("The measurement took {}", elapsed);
    }

    TEST_CASE("Vec3d: operations - project / reject (vector - bivector)")
    {
        fmt::println("Vec3d: operations - project / reject (vector - bivector)");

        Vec3d v1{5.0, 3.0, 1.0};
        BiVec3d v2 = wdg(Vec3d{0.0, 0.0, 2.0}, Vec3d{2.0, 0.0, 0.0});

        Vec3d v3{project_onto(v1, v2)};
        Vec3d v4{reject_from(v1, v2)};
        Vec3d v5{v3 + v4};

        // fmt::println("v1  = {: .8f}, nrm(v1) = {: .8f}", v1, nrm(v1));
        // fmt::println("v2  = {: .8f}, nrm(v2) = {: .8f}", v2, nrm(v2));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .8f}, nrm(v3) = {: .8f}", v3,
        //              nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .8f}, nrm(v4) = {: .8f}", v4,
        //              nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .8f}, nrm(v5) = {: .8f}", v5,
        //              nrm(v5));

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);
    }

    TEST_CASE("Vec3d: cross-product")
    {
        fmt::println("Vec3d: cross-product");

        Vec3d u{1.0, 1.0, 0.0};
        Vec3d v{0.0, 1.0, 1.0};
        Vec3d w{1.0, 1.0, 1.0};

        Vec3d u_cross_v = cross(u, v);
        BiVec3d u_wdg_v = wdg(u, v);
        auto minus_dual_u_wdg_v = -dual(wdg(u, v));

        // fmt::println("   u          = {: } - vec,   nrm(u)"
        //              "          = {: .4}",
        //              u, nrm(u));
        // fmt::println("   v          = {: } - vec,   nrm(v)"
        //              "          = {: .4}",
        //              v, nrm(v));
        // fmt::println("   uxv        = {: } - vec,   nrm(uxv)"
        //              "        = {: .4}",
        //              u_cross_v, nrm(u_cross_v));
        // fmt::println("   u^v        = {: } - bivec, nrm(u^v)"
        //              "        = {: .4}",
        //              u_wdg_v, nrm(u_wdg_v));
        // fmt::println("   -dual(u^v) = {: } - vec  , nrm(-dual(u^v))"
        //              " = {: .4}",
        //              minus_dual_u_wdg_v, nrm(minus_dual_u_wdg_v));
        // fmt::println("");
        // fmt::println("   u                = {: } - vec, nrm(u)"
        //              "                = {: .4}",
        //              u, nrm(u));
        // fmt::println("   v                = {: } - vec, nrm(v)"
        //              "                = {: .4}",
        //              v, nrm(v));
        // fmt::println("   w                = {: } - vec, nrm(w)"
        //              "                = {: .4}",
        //              w, nrm(w));
        // fmt::println("   u x (v x w)      = {: } - vec, nrm(u x v x w)"
        //              "        = {: .4}",
        //              cross(u, cross(v, w)), nrm(cross(u, cross(v, w))));
        // fmt::println("   -dot(u, wdg(v,w) = {: } - vec, nrm(-dot(u, wdg(v,w))"
        //              " = {: .4}",
        //              -dot(u, wdg(v, w)), nrm(-dot(u, wdg(v, w))));
        // fmt::println("");

        CHECK(u_cross_v == minus_dual_u_wdg_v);
        CHECK(cross(u, cross(v, w)) == -dot(u, wdg(v, w)));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d: default init")
    {
        fmt::println("MVec3d: default init");
        // default initialization
        MVec3d v;
        // fmt::println("   v = {}", v);
        CHECK(v.c0 == 0.0);
        CHECK(v.c1 == 0.0);
        CHECK(v.c2 == 0.0);
        CHECK(v.c3 == 0.0);
        CHECK(v.c4 == 0.0);
        CHECK(v.c5 == 0.0);
        CHECK(v.c6 == 0.0);
        CHECK(v.c7 == 0.0);
    }
    TEST_CASE("MVec3d: with curly braced intializer")
    {
        fmt::println("MVec3d: with curly braced intializer");
        // default initialization
        MVec3d v{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        // fmt::println("   v = {}", v);
        CHECK(v.c0 == 0.0);
        CHECK(v.c1 == 1.0);
        CHECK(v.c2 == 2.0);
        CHECK(v.c3 == 3.0);
        CHECK(v.c4 == 23.0);
        CHECK(v.c5 == 31.0);
        CHECK(v.c6 == 12.0);
        CHECK(v.c7 == 123.0);
    }

    TEST_CASE("MVec3d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("MVec3d: cp ctor & cp assign incl. type deduction");
        // default initialization
        MVec3d v1{0.0,  1.0,  2.0,  3.0,
                  23.0, 31.0, 12.0, 123.0}; // init with double (type deduction)
        MVec3d v2{v1};                      // cp ctor
        MVec3d v3 = v2;                     // cp assign
        MVec3d v4 = -v3;                    // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(v2.c0 == 0.0);
        CHECK(v2.c1 == 1.0);
        CHECK(v2.c2 == 2.0);
        CHECK(v2.c3 == 3.0);
        CHECK(v2.c4 == 23.0);
        CHECK(v2.c5 == 31.0);
        CHECK(v2.c6 == 12.0);
        CHECK(v2.c7 == 123.0);

        CHECK(v3.c0 == 0.0);
        CHECK(v3.c1 == 1.0);
        CHECK(v3.c2 == 2.0);
        CHECK(v3.c3 == 3.0);
        CHECK(v3.c4 == 23.0);
        CHECK(v3.c5 == 31.0);
        CHECK(v3.c6 == 12.0);
        CHECK(v3.c7 == 123.0);
    }

    TEST_CASE("MVec3d: fmt & cout printing")
    {
        fmt::println("MVec3d: fmt & cout printing");

        MVec3d pf{1.0f, 2.0001f, 0.0f, 3.0f, 1.0f, 2.0001f, 0.0f, 3.0f};
        MVec3d pd{1.0, 2.0001, 0.0, 3.0, 1.0, 2.0001, 0.0, 3.0};

        // std::cout << "    cout: pf = " << pf << std::endl;
        // fmt::println("    fmt:  pf = {}", pf);
        // fmt::println("    fmt:  pf = {:.8f}", pf);

        // std::cout << "    cout: pd = " << pd << std::endl;
        // fmt::println("    fmt:  pd = {}", pd);
        // fmt::println("    fmt:  pd = {:.8f}", pd);

        // std::vector<MVec3d<double>> vp1{{1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0},
        //                                 {0.5, 1.5, 2.0, 2.5, 1.0, 1.0, 1.0, 2.0}};
        // fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        // fmt::println("    fmt: vp1 = {:.e}", fmt::join(vp1, ", "));

        CHECK(pf == pd);
    }

    TEST_CASE("MVec3d: vector space and linearity tests")
    {
        fmt::println("MVec3d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        MVec3d p0;
        MVec3d p1{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
        MVec3d p2{0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0};
        MVec3d p3{0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0};
        MVec3d p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t.multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t. addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d: geometric product tests - gpr(vec,vec)")
    {
        fmt::println("MVec3d: geometric product tests - gpr(vec,vec)");

        // ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
        //
        // dot(a,b) = 0.5*(ab + ba)   (symmetric part)
        // wdg(a,b) = 0.5*(ab - ba)   (antisymmetric part)

        Vec3d a{1.0, 2.0, 3.0};
        Vec3d b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        MVec3d mva{a};
        MVec3d mvb{b};
        auto mvab = gpr(mva, mvb);
        auto mvab_sym = 0.5 * (gpr(mva, mvb) + gpr(mvb, mva));
        auto mvab_asym = 0.5 * (gpr(mva, mvb) - gpr(mvb, mva));

        // fmt::println("   a = {}", a);
        // fmt::println("   b = {}", b);
        // fmt::println("   dot(a,b) = {}", dot_ab);
        // fmt::println("   wdg(a,b) = {}", wdg_ab);
        // fmt::println("");
        // fmt::println("   mva = {}", mva);
        // fmt::println("   mvb = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym = 0.5*(mva mvb + mvb mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva mvb - mvb mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));

        CHECK(dot_ab == gr0(mvab));
        CHECK(dot_ab == gr0(mvab_sym));
        CHECK(wdg_ab == gr2(mvab));
        CHECK(wdg_ab == gr2(mvab_asym));
    }

    TEST_CASE("MVec3d: geometric product tests - gpr(bivec,vec)")
    {
        fmt::println("MVec3d: geometric product tests - gpr(bivec,vec)");

        // Ab = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)
        //
        // dot(A,b) = 0.5*(Ab - Aa)   (antisymmetric part)
        // wdg(A,b) = 0.5*(Ab + Aa)   (symmetric part)

        BiVec3d a{1.0, 2.0, 3.0};
        Vec3d b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        MVec3d mva{a};
        MVec3d mvb{b};
        auto mvab = gpr(mva, mvb);
        auto mvab_sym = 0.5 * (gpr(mva, mvb) + gpr(mvb, mva));
        auto mvab_asym = 0.5 * (gpr(mva, mvb) - gpr(mvb, mva));

        // fmt::println("   a = {}", a);
        // fmt::println("   b = {}", b);
        // fmt::println("   dot(a,b) = {}", dot_ab);
        // fmt::println("   wdg(a,b) = {}", wdg_ab);
        // fmt::println("");
        // fmt::println("   mva = {}", mva);
        // fmt::println("   mvb = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym = 0.5*(mva mvb + mvb mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva mvb - mvb mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));

        CHECK(dot_ab == gr1(mvab));
        CHECK(dot_ab == gr1(mvab_asym));
        CHECK(wdg_ab == gr3(mvab));
        CHECK(wdg_ab == gr3(mvab_sym));
    }

    TEST_CASE("MVec3d: geometric product tests - gpr(vec,bivec)")
    {
        fmt::println("MVec3d: geometric product tests - gpr(vec,bivec)");

        // aB = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)
        //
        // dot(a,B) = 0.5*(aB - Ba)   (antisymmetric part)
        // wdg(a,B) = 0.5*(aB + Ba)   (symmetric part)

        Vec3d a{1.0, 2.0, 3.0};
        BiVec3d b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        MVec3d mva{a};
        MVec3d mvb{b};
        auto mvab = gpr(mva, mvb);
        auto mvab_sym = 0.5 * (gpr(mva, mvb) + gpr(mvb, mva));
        auto mvab_asym = 0.5 * (gpr(mva, mvb) - gpr(mvb, mva));

        // fmt::println("   a = {}", a);
        // fmt::println("   b = {}", b);
        // fmt::println("   dot(a,b) = {}", dot_ab);
        // fmt::println("   wdg(a,b) = {}", wdg_ab);
        // fmt::println("");
        // fmt::println("   mva = {}", mva);
        // fmt::println("   mvb = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym = 0.5*(mva mvb + mvb mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva mvb - mvb mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));

        CHECK(dot_ab == gr1(mvab));
        CHECK(dot_ab == gr1(mvab_asym));
        CHECK(wdg_ab == gr3(mvab));
        CHECK(wdg_ab == gr3(mvab_sym));
    }

    TEST_CASE("MVec3d geometric product tests - recovering vectors from the"
              "geometric product")
    {
        fmt::println("MVec3d: geometric product tests - recovering vectors from the "
                     "geometric product");

        // Two multivectors mv1 and mv2 formed from vectors v1 and v2.
        // (gr0(mv1)==0 && gr1(mv1) != 0 && gr2(mv1)==0 &&
        //  gr0(mv2)==0 && gr1(mv2) != 0 && gr2(mv2)==0 )
        //
        // They are multiplied by the geometric product to form a multivector C
        // C = mv1(v1) mv2(v2) = gpr(mv1,mv2)
        //
        // C contains a scalar part and a bivector part exclusively,
        // the remaining components are zero.
        // (gr0(C) != 0 && gr1(C)==0 && gr2(C) !=0)
        //
        // The scalar part of C represents the parts of v1 and v2
        // that are parallel to each other.
        // The bivector part of C represents the parts of v1 and v2
        // that are perpendicular to each other.
        //
        // multiply C from the right with inv(v2) recovers v1
        // multiply C from the left the the inv(v1) recovers v2

        Vec3d a{1.0, 2.0, 3.0};
        Vec3d b{0.5, 3.0, -4.0};
        MVec3d mva{a};
        MVec3d mvb{b};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);
        MVec3d C = gpr(a, b);
        MVec3d Cm = gpr(mva, mvb);
        MVec3d Cd{Scalar3d<double>(dot_ab), wdg_ab};

        MVec3d gpr_right = gpr(C, MVec3d{inv(b)});
        MVec3d gpr_left = gpr(MVec3d{inv(a)}, C);

        // fmt::println("   a                           = {}", a);
        // fmt::println("   b                           = {}", b);
        // fmt::println("   C  = gpr(a,b)               = {}", C);
        // fmt::println("   Cm = gpr(mva,mvb)           = {}", Cm);
        // fmt::println("   Cd = mv(dot(a,b), wdg(a,b)) = {}", Cd);
        // fmt::println("");
        // fmt::println("   gpr(C,bi) = gpr_right = {}", gpr_right);
        // fmt::println("   gpr(ai,C) = gpr_left = {}", gpr_left);
        // fmt::println("   gr1(gpr_right) = a = {}", gr1(gpr_right));
        // fmt::println("   gr1(gpr_left) = b = {}", gr1(gpr_left));

        CHECK(C == Cm);
        CHECK(C == Cd);
        CHECK(a == gr1(gpr_right));
        CHECK(b == gr1(gpr_left));
    }

    TEST_CASE("MVec3d: geometric product tests - equivalence tests")
    {
        fmt::println("MVec3d: geometric product tests - equivalence tests");

        Vec3d a{1.0, 2.0, 3.0};
        Vec3d b{0.5, 3.0, -4.0};
        MVec3d mva{a};
        MVec3d mvb{b};

        BiVec3d A{1.0, 2.0, 3.0};
        BiVec3d B{0.5, 3.0, -4.0};
        MVec3d mvA{A};
        MVec3d mvB{B};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        auto dot_Ab = dot(A, b);
        auto wdg_Ab = wdg(A, b);

        auto dot_aB = dot(a, B);
        auto wdg_aB = wdg(a, B);

        MVec3d_E ab = gpr(a, b);
        MVec3d abm = gpr(mva, mvb);
        MVec3d abd{Scalar3d<double>(dot_ab), wdg_ab};

        MVec3d_U Ab = gpr(A, b);
        MVec3d Abm = gpr(mvA, mvb);
        MVec3d Abd{dot_Ab, wdg_Ab};

        MVec3d_U aB = gpr(a, B);
        MVec3d aBm = gpr(mva, mvB);
        MVec3d aBd{dot_aB, wdg_aB};

        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   ab  = MVec3d_E(gpr(a,b))         = {}", ab);
        // fmt::println("   abm = gpr(mva,mvb)               = {}", abm);
        // fmt::println("   abd = MVec3d(dot(a,b), wdg(a,b)) = {}", abd);
        // fmt::println("");
        // fmt::println("   A                                = {}", A);
        // fmt::println("   mvA                              = {}", mvA);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   Ab  = MVec3d_U(gpr(A,b))         = {}", Ab);
        // fmt::println("   Abm = gpr(mvA,mvb)               = {}", Abm);
        // fmt::println("   Abd = MVec3d(dot(A,b), wdg(A,b)) = {}", Abd);
        // fmt::println("");
        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   B                                = {}", B);
        // fmt::println("   mvB                              = {}", mvB);
        // fmt::println("   aB  = MVec3d_U(gpr(a,B))         = {}", aB);
        // fmt::println("   aBm = gpr(mva,mvB)               = {}", aBm);
        // fmt::println("   aBd = MVec3d(dot(a,B), wdg(a,B)) = {}", aBd);
        // fmt::println("");

        CHECK(gr0(ab) == gr0(abm));
        CHECK(gr1(abm) == Vec3d{});
        CHECK(gr2(ab) == gr2(abm));
        CHECK(gr3(abm) == PScalar3d<double>{0.0});

        CHECK(gr0(ab) == gr0(abd));
        CHECK(gr1(abd) == Vec3d{});
        CHECK(gr2(ab) == gr2(abd));
        CHECK(gr3(abd) == PScalar3d<double>{0.0});

        CHECK(gr0(Abm) == 0);
        CHECK(gr1(Ab) == gr1(Abm));
        CHECK(gr2(Abm) == BiVec3d{});
        CHECK(gr3(Ab) == gr3(Abm));

        CHECK(gr0(Abd) == 0);
        CHECK(gr1(Ab) == gr1(Abd));
        CHECK(gr2(Abd) == BiVec3d{});
        CHECK(gr3(Ab) == gr3(Abd));

        CHECK(gr0(aBm) == 0);
        CHECK(gr1(aB) == gr1(aBm));
        CHECK(gr2(aBm) == BiVec3d{});
        CHECK(gr3(aB) == gr3(aBm));

        CHECK(gr0(aBd) == 0);
        CHECK(gr1(aB) == gr1(aBd));
        CHECK(gr2(aBd) == BiVec3d{});
        CHECK(gr3(aB) == gr3(aBd));
    }

    TEST_CASE("MVec3d: assignment tests")
    {
        fmt::println("MVec3d: assignment tests");

        Vec3d v1{1.0, 2.0, 3.0};
        Vec3d v2{0.5, 1.0, 1.5};
        Vec3d v3{0.5, 1.0, -4.5};
        BiVec3d b1{1.0, 2.0, 3.0};

        MVec3d<double> mv0;
        MVec3d mv1{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        MVec3d mv2{0.0, 0.5, 1.0, 1.5, 11.5, 15.5, 6.0, 61.5};
        MVec3d mv3{mv1};
        MVec3d mv4 = mv2;

        MVec3d mv5(Scalar3d<double>(5.0));
        MVec3d mv6(PScalar3d<double>(6.0));
        MVec3d mv7{v1};
        MVec3d mv8{b1};
        MVec3d mv9{Scalar3d<double>(dot(v1, v3)), wdg(v1, v3)};

        MVec3d mv10{v1, PScalar3d<double>(10.0)};
        // This must not compile! Implict conversion to Vec3d possible
        // possible solution: explicitly deleted constructor for MVec3d
        // MVec3d mv11{b1, pscalar3d_t(10.0)};

        // this does not compile (which is fine, a base cannot convert to derived)
        // MVec3d mv12{Scalar3d<double>(10.0), v1};

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("");
        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("   mv3 = {}", mv3);
        // fmt::println("   mv4 = {}", mv4);
        // fmt::println("   mv5 = {}", mv5);
        // fmt::println("   mv6 = {}", mv6);
        // fmt::println("");
        // fmt::println("   gr1(mv1) = {}", gr1(mv1));
        // fmt::println("   gr1(mv2) = {}", gr1(mv2));
        // fmt::println("   gr1(mv3) = {}", gr1(mv3));
        // fmt::println("   gr1(mv3) = {}", gr1(mv4));
        // fmt::println("");
        // fmt::println("   v1 = {}", v1);
        // fmt::println("   mv7 = v1 = {}", mv7);
        // fmt::println("   b1 = {}", b1);
        // fmt::println("   mv8 = b1 = {}", mv8);
        // fmt::println("");
        // fmt::println("   mv9 = {}", mv9);
        // fmt::println("   mv10 = {}", mv10);

        CHECK(gr1(mv1) == v1);
        CHECK(gr1(mv2) == v2);
        CHECK(gr1(mv3) == v1);
        CHECK(gr1(mv4) == v2);
        CHECK(gr0(mv5) == 5.0);
        CHECK(gr3(mv6) == 6.0);
        CHECK(mv1 == mv3);
        CHECK(mv4 == mv2);
        CHECK(gr1(mv7) == v1);
        CHECK(gr2(mv8) == b1);
        CHECK(gr0(mv9) == dot(v1, v3));
        CHECK(gr2(mv9) == wdg(v1, v3));
        CHECK(gr1(mv10) == v1);
        CHECK(gr3(mv10) == 10.0);
    }

    TEST_CASE("MVec3d: bivector product properties")
    {
        fmt::println("MVec3d: bivector product properties");

        BiVec3d b1{1.0, 2.0, 3.0};
        MVec3d mb1{b1};
        BiVec3d b2{-3.0, 1.0, 2.0};
        MVec3d mb2{b2};

        auto b12_dot = dot(b1, b2);
        auto b12_cmt = cmt(b1, b2);
        auto b21_dot = dot(b2, b1);
        auto b21_cmt = cmt(b2, b1);

        auto gpr12_m = gpr(mb1, mb2);
        auto gpr21_m = gpr(mb2, mb1);
        auto gpr12_m_sym = 0.5 * (gpr12_m + gpr21_m);
        auto gpr12_m_asym = 0.5 * (gpr12_m - gpr21_m);

        auto gpr12_d = gpr(b1, b2);
        auto gpr21_d = gpr(b2, b1);
        auto gpr12_d_sym = 0.5 * (gpr12_d + gpr21_d);
        auto gpr12_d_asym = 0.5 * (gpr12_d - gpr21_d);

        // fmt::println("   b1  = {}", b1);
        // fmt::println("   mb1 = {}", mb1);
        // fmt::println("   b2  = {}", b2);
        // fmt::println("   mb2 = {}", mb2);
        // fmt::println("");
        // fmt::println("   b12_dot = {}", b12_dot);
        // fmt::println("   b12_cmt = {}", b12_cmt);
        // fmt::println("   b21_dot = {}", b21_dot);
        // fmt::println("   b21_cmt = {}", b21_cmt);
        // fmt::println("");
        // fmt::println("   gpr12_m = gpr(mb1, mb2) = {}", gpr12_m);
        // fmt::println("   gpr21_m = gpr(mb2, mb1) = {}", gpr21_m);
        // fmt::println("   gpr12_m_sym  = 0.5 * (gpr12_d + gpr21_d) = {}",
        // gpr12_m_sym); fmt::println("   gpr12_m_asym = 0.5 * (gpr12_m - gpr21_m) =
        // {}", gpr12_m_asym); fmt::println(""); fmt::println("   gpr12_d = gpr(b1,
        // b2) = {}", gpr12_d); fmt::println("   gpr21_d = gpr(b2, b1) = {}",
        // gpr21_d); fmt::println("   gpr12_d_sym  = 0.5 * (gpr12_d + gpr21_d) = {}",
        // gpr12_d_sym); fmt::println("   gpr12_d_asym = 0.5 * (gpr12_d - gpr21_d) =
        // {}", gpr12_d_asym); fmt::println("");

        CHECK(gr2(mb1) == b1);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3d_E<T> and MVec3d_U<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d_E/_U: modelling even and uneven parts of 3d algebra - basics")
    {
        fmt::println(
            "MVec3d_E/_U: modelling even and uneven parts of 3d algebra - basics");

        // defining a complex number in all three forms as multivector
        auto u = unitized(Vec3d{1.0, 0.0, 0.0});
        auto v =
            unitized(Vec3d{std::cos(pi / 12), std::sin(pi / 12), 0.0}); // unit vec +15%
        auto angle_uv = angle(u, v);
        auto B = wdg(u, v); // unitized bivector describing the plane spanned by u and v

        auto my_exp = exp(-B, angle_uv);
        auto my_rot = rotor(B, 2.0 * angle_uv);

        // definition of rotor used here: B = u^v
        // => B determines the meaning of the positive sign of the rotation
        //
        auto R_m =
            MVec3d(exp(-B, angle_uv)); // Rotor formed by u and v (unitized bivector)
        auto Rr_m = MVec3d(rev(R_m));  // and its reverse

        auto c = Vec3d{1.0, 1.0, 1.0};
        auto c_m = MVec3d{c};

        auto c_tmp_m = R_m * c_m;
        auto c_rot_m = c_tmp_m * Rr_m;

        auto R = exp(-B, angle_uv); // Rotor formed by u and v (unitized bivector)
        auto Rr = rev(R);           // and its reverse

        auto c_tmp_l = R * c;
        auto c_rot_u_l = c_tmp_l * Rr;
        auto c_rot_l = gr1(c_rot_u_l);
        // due to symmetry of R and Rr the gr3(c_rot) part will be zero
        // and thus can be assumed to be zero for further computations

        auto c_tmp_r = c * Rr;
        auto c_rot_u_r = R * c_tmp_r;
        auto c_rot_r = gr1(c_rot_u_r);
        // due to symmetry of R and Rr the gr3(c_rot) part will be zero
        // and thus can be assumed to be zero for further computations

        auto angle_c_c_rot = angle(c, c_rot_l); // not that easy in 3D!
        // (angle in plane of both vectors is not the angle in the plane
        // represented by the bivector!)
        // => requires projection of vectors onto plane and then taking
        // the angle between the projected vectors to be correct (bivector angle!)

        auto c_proj = project_onto(c, B);
        auto c_rot_proj = project_onto(c_rot_l, B);
        auto angle_proj = angle(c_proj, c_rot_proj);

        // fmt::println("   u                     = {: .3}", u);
        // fmt::println("   v                     = {: .3}", v);
        // fmt::println("   B = u^v = wdg(u,v)    = {: .3}", B);
        // fmt::println("   angle(u,v)            = {: .3}°", rad_to_deg(angle_uv));
        // fmt::println("   sin(angle(u,v))       = {: .3}", std::sin(angle_uv));
        // fmt::println("");
        // fmt::println("   c                     = {: .3}", c);
        // fmt::println("");
        // fmt::println("Implemented as full multivector operation:");
        // fmt::println("   R_m  = MVec3d(exp(-B,angle_uv))  = {: .3}", R_m);
        // fmt::println("   Rr_m = rev(R_m)                  = {: .3}", Rr_m);
        // fmt::println("   Rr_m*R_m                         = {: .3}", Rr_m * R_m);
        // fmt::println("   c_m                              = {: .3}", c_m);
        // fmt::println("   c_tmp_m = R_m*c_m                = {: .3}", c_tmp_m);
        // fmt::println("   c_rot_m = c_tmp_m*Rr_m           = {: .3}", c_rot_m);
        // fmt::println("   gr1(c_rot_m)                     = {: .3}", gr1(c_rot_m));
        // fmt::println("");
        // fmt::println("Implemented as reduced grade multivector operation:");
        // fmt::println("   R  = exp(-B,angle_uv)            = {: .3}", R);
        // fmt::println("   Rr = rev(R)                      = {: .3}", Rr);
        // fmt::println("   my_exp = exp(-B, angle_uv)       = {: .3}", my_exp);
        // fmt::println("   my_rot = rotor(B, 2*angle_uv)    = {: .3}", my_rot);
        // fmt::println("");
        // fmt::println("Left multiplication of rotor first:");
        // fmt::println("   c_tmp_l = R*c            = {: .3}", c_tmp_l);
        // fmt::println("   c_rot_u_l = c_tmp_l*Rr   = {: .3}", c_rot_u_l);
        // fmt::println("   c_rot_l = gr1(c_rot_u_l) = {: .3}", c_rot_l);
        // fmt::println("");
        // fmt::println("Right multiplication of rotor first:");
        // fmt::println("   c_tmp_r = c*Rr           = {: .3}", c_tmp_r);
        // fmt::println("   c_rot_u_r = R*c_tmp_r    = {: .3}", c_rot_u_r);
        // fmt::println("   c_rot_r = gr1(c_rot_u_r) = {: .3}", c_rot_r);
        // fmt::println("");
        // fmt::println("   angle(c, c_rot_l) = {: .3}°", rad_to_deg(angle_c_c_rot));
        // fmt::println("   angle(projected)  = {: .3}°", rad_to_deg(angle_proj));
        // fmt::println("");
        // fmt::println("direct calclulation:");
        // fmt::println("   c_rot = rotate(c,R)          = {: .3}", rotate(c, R));

        CHECK(nrm(rotate(c, R)) == nrm(c));
        CHECK(gr1(c_rot_m) == rotate(c, R));
        // n I_3d approach:
        CHECK(rotate(Vec3d{1.0, 0.0, 0.0}, rotor(e3_3d * I_3d, pi / 4)) ==
              unitized(Vec3d{1.0, 1.0, 0.0}));
        // using a bivector directly:
        CHECK(rotate(Vec3d{1.0, 0.0, 0.0}, rotor(e12_3d, pi / 4)) ==
              unitized(Vec3d{1.0, 1.0, 0.0}));

        // direct rotation of a bivector
        CHECK(rotate(BiVec3d{0.0, 0.0, 1.0}, rotor(e23_3d, pi / 2)) == -e31_3d);
    }

    TEST_CASE("MVec3d: dualization")
    {
        fmt::println("MVec3d: dualization");


        Vec3d v{1.0, 2.0, 3.0};                                   // 3d vector
        BiVec3d B{10.0, 20.0, 30.0};                              // 3d bivector
        MVec3d vm{100.0, 1.0, 2.0, 3.0, 10.0, 20.0, 30.0, 300.0}; // full 3d multivector

        // full 3d multivector - even content
        MVec3d vm_even{100.0, 0.0, 0.0, 0.0, 10.0, 20.0, 30.0, 0.0};
        // even grade 3d multivector
        MVec3d_E vm_E{100.0, 10.0, 20.0, 30.0};

        // full 3d multivector - uneven content
        MVec3d vm_uneven{0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 300.0};
        // uneven grade 3d multivector
        MVec3d_U vm_U{1.0, 2.0, 3.0, 300.0};

        // dual by left multiplication with Im_3d
        // as defined in Doran/Lasenby "GA for physicists"

        auto vm_dual_manual = Im_3d * vm;
        auto vm_dual = dual(vm);

        auto vm_dual_even_manual = Im_3d * vm_even;
        auto vm_dual_even = dual(vm_even);

        auto vm_dual_uneven_manual = Im_3d * vm_uneven;
        auto vm_dual_uneven = dual(vm_uneven);

        // result is uneven, naming chosen for consistency
        auto vm_dual_manual_E = I_3d * vm_E;
        auto vm_dual_E = dual(vm_E);

        // result is even, naming chosen for consistency
        auto vm_dual_manual_U = Im_3d_U * vm_U;
        auto vm_dual_U = dual(vm_U);

        auto v_dual_manual = I_3d * v;
        auto v_dual = dual(v);

        auto B_dual_manual = I_3d * B;
        auto B_dual = dual(B);

        // fmt::println("   I_3d    = {}", I_3d);
        // fmt::println("   Im_3d   = {}", Im_3d);
        // fmt::println("   Im_3d_U = {}", Im_3d_U);
        // fmt::println("");
        // fmt::println("   v             = {}", v);
        // fmt::println("   B             = {}", B);
        // fmt::println("");
        // fmt::println("   vm            = {}", vm);
        // fmt::println("   Im_3d*vm      = {}", vm_dual_manual);
        // fmt::println("   dual(vm)      = {}", vm_dual);
        // fmt::println("");
        // fmt::println("   vm_even       = {}", vm_even);
        // fmt::println("   Im_3d*vm_even = {}", vm_dual_even_manual);
        // fmt::println("   dual(vm_even) = {}", vm_dual_even);
        // fmt::println("");
        // fmt::println("   vm_E          = {}", vm_E);
        // fmt::println("   Im_3d_E*vm_E  = {}", vm_dual_manual_E);
        // fmt::println("   dual(vm_E)    = {}", vm_dual_E);
        // fmt::println("");
        // fmt::println("   vm_uneven       = {}", vm_uneven);
        // fmt::println("   Im_3d*vm_uneven = {}", vm_dual_uneven_manual);
        // fmt::println("   dual(vm_uneven) = {}", vm_dual_uneven);
        // fmt::println("");
        // fmt::println("   vm_U          = {}", vm_U);
        // fmt::println("   Im_3d_U*vm_U  = {}", vm_dual_manual_U);
        // fmt::println("   dual(vm_U)    = {}", vm_dual_U);
        // fmt::println("");
        // fmt::println("   v               = {}", v);
        // fmt::println("   I_3d * v        = {} - bivec ", v_dual_manual);
        // fmt::println("   dual(v)         = {} - bivec ", v_dual);
        // fmt::println("");
        // fmt::println("   B               = {}", B);
        // fmt::println("   I_3d * B        = {} - vec", B_dual_manual);
        // fmt::println("   dual(B)         = {} - vec", B_dual);

        CHECK(vm_dual == vm_dual_manual);
        CHECK(vm_dual_even == vm_dual_even_manual);
        CHECK(vm_dual_uneven == vm_dual_uneven_manual);
        CHECK(vm_dual_E == vm_dual_manual_E);
        CHECK(vm_dual_U == vm_dual_manual_U);
        CHECK(dual(v) == BiVec3d{1.0, 2.0, 3.0});
        CHECK(dual(B) == -Vec3d{10.0, 20.0, 30.0});
        CHECK(dual(Scalar3d<double>(5)) == PScalar3d<double>(5));
        CHECK(dual(PScalar3d<double>(5)) == Scalar3d<double>(-5));
    }

} // TEST_SUITE("Geometric Algebra")
