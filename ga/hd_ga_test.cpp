#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <chrono>
#include <numbers> // math constants like pi_v
#include <vector>

// include functions to be tests
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

    TEST_CASE("Vec2d default init")
    {
        fmt::println("");
        fmt::println("Vec2d default init:");
        // default initialization
        Vec2d v;
        CHECK(v.x == 0.0f);
        CHECK(v.y == 0.0f);
        fmt::println("   v = {}", v);
    }
    TEST_CASE("Vec2d with curly braced intializer")
    {
        fmt::println("");
        fmt::println("Vec2d with curly braced intializer:");
        // default initialization
        Vec2d v{0.0f, 0.0f};
        CHECK(v.x == 0.0f);
        CHECK(v.y == 0.0f);
        fmt::println("   v = {}", v);
    }
    TEST_CASE("Vec2d cp ctor & cp assign incl. type deduction")
    {
        fmt::println("");
        fmt::println("Vec2d cp ctor & cp assign incl. type deduction:");
        // default initialization
        Vec2d v1{1.0, 2.0}; // init with double (type deduction)
        Vec2d v2{v1};       // cp ctor
        Vec2d v3 = v2;      // cp assign
        CHECK(v1.x == 1.0);
        CHECK(v1.y == 2.0);
        CHECK(v2.x == 1.0);
        CHECK(v2.y == 2.0);
        CHECK(v3.x == 1.0);
        CHECK(v3.y == 2.0);
        fmt::println("   v1 = {}", v1);
        fmt::println("   v2 = {}", v2);
        fmt::println("   v3 = {}", v3);
    }


    TEST_CASE("Vec2d fmt & cout printing")
    {
        fmt::println("");
        fmt::println("Vec2d fmt & cout printing:");

        Vec2d pf{1.0f, 2.0001f};
        Vec2d pd{1.0, 2.0001};

        std::cout << "   cout: pf = " << pf << std::endl;
        fmt::println("    fmt: pf = {}", pf);
        fmt::println("    fmt: pf = {:.7f}", pf);

        std::cout << "   cout: pd = " << pd << std::endl;
        fmt::println("    fmt: pd = {}", pd);
        fmt::println("    fmt: pd = {:.7f}", pd);

        std::vector<Vec2d<double>> vp1{{1.0, 1.0}, {1.5, 2.0}};
        fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("    fmt: vp1 = {:.e}", fmt::join(vp1, ", "));
    }

    TEST_CASE("Vec2d comparison float")
    {
        fmt::println("");
        fmt::println("Vec2d comparison float:");

        Vec2d<float> v1f{1.0f, 2.0f};
        Vec2d<float> v2f{2.0f, 4.0f};
        Vec2d<float> v3f{1.0f, 2.0000001f};
        Vec2d<float> v4f{v1f};

        fmt::println("   v1f = {}", v1f);
        fmt::println("   v2f = {}", v2f);
        fmt::println("   v3f = {}", v3f);
        fmt::println("   v4f = {}", v4f);

        fmt::println("    fmt: eps = {}", std::numeric_limits<float>::epsilon());

        CHECK(v1f == v4f);             // comparison (equality)
        CHECK(v1f != v2f);             // comparison (inequality)
        CHECK(norm(v1f) < norm(v2f));  // comparison (less than)
        CHECK(norm(v2f) >= norm(v1f)); // comparison (greater than or equal)
        CHECK(v3f == v1f);             // comparison (eqality)
    }

    TEST_CASE("Vec2d comparison double")
    {
        fmt::println("");
        fmt::println("Vec2d comparison double:");

        Vec2d<double> v1d{1.0, 2.0};
        Vec2d<double> v2d{2.0, 4.0};
        Vec2d<double> v3d{1.0, 2.0000000000000001};
        Vec2d<double> v4d{v1d};

        fmt::println("   v1d = {}", v1d);
        fmt::println("   v2d = {}", v2d);
        fmt::println("   v3d = {}", v3d);
        fmt::println("   v4d = {}", v4d);

        fmt::println("    fmt: eps = {}", std::numeric_limits<double>::epsilon());

        CHECK(v1d == v4d);             // comparison (equality)
        CHECK(v1d != v2d);             // comparison (inequality)
        CHECK(norm(v1d) < norm(v2d));  // comparison norm
        CHECK(norm(v2d) >= norm(v1d)); // comparison norm
        CHECK(v3d == v1d);             // comparison (eqality)
    }

    TEST_CASE("Vec2d operations - norm, inverse, dot")
    {
        fmt::println("");
        fmt::println("Vec2d operations - norm, inverse, dot:");

        vec2d v1{2.0, 1.0};
        vec2d v2{normalized(v1)};

        vec2d v3{2.0, 6.0};
        vec2d v4{inverse(v3)};

        fmt::println("v1 = {: .5f}, norm(v1) = {: .5f}", v1, norm(v1));
        fmt::println("v2 = normalized(v1) = {: .5f}, norm(v2) = {: .5f}", v2, norm(v2));

        CHECK(std::abs(sq_norm(v1) - 5.0) < eps);
        CHECK(std::abs(sq_norm(v2) - 1.0) < eps);
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);
    }

    TEST_CASE("Vec2d operations - angle")
    {
        fmt::println("");
        fmt::println("Vec2d operations - angle:");

        vec2d v1{1.0, 0.0};
        vec2d v2{normalized(vec2d(1.0, 1.0))};
        vec2d v3{0.0, 1.0};
        vec2d v4{normalized(vec2d(-1.0, 1.0))};
        vec2d v5{-1.0, 0.0};
        vec2d v6{normalized(vec2d(-1.0, -1.0))};
        vec2d v7{0.0, -1.0};
        vec2d v8{normalized(vec2d(1.0, -1.0))};

        using std::numbers::pi;

        fmt::println("v1 = {: .5f}, norm(v1) = {: .5f}, angle(v1,v1) = {: .5f}, {: .5f}",
                     v1, norm(v1), angle(v1, v1), angle(v1, v1) / pi);
        fmt::println("v2 = {: .5f}, norm(v2) = {: .5f}, angle(v1,v2) = {: .5f}, {: .5f}",
                     v2, norm(v2), angle(v1, v2), angle(v1, v2) / pi);
        fmt::println("v3 = {: .5f}, norm(v3) = {: .5f}, angle(v1,v3) = {: .5f}, {: .5f}",
                     v3, norm(v3), angle(v1, v3), angle(v1, v3) / pi);
        fmt::println("v4 = {: .5f}, norm(v4) = {: .5f}, angle(v1,v4) = {: .5f}, {: .5f}",
                     v4, norm(v4), angle(v1, v4), angle(v1, v4) / pi);
        fmt::println("v5 = {: .5f}, norm(v5) = {: .5f}, angle(v1,v5) = {: .5f}, {: .5f}",
                     v5, norm(v5), angle(v1, v5), angle(v1, v5) / pi);
        fmt::println("v6 = {: .5f}, norm(v6) = {: .5f}, angle(v1,v6) = {: .5f}, {: .5f}",
                     v6, norm(v6), angle(v1, v6), angle(v1, v6) / pi);
        fmt::println("v7 = {: .5f}, norm(v7) = {: .5f}, angle(v1,v7) = {: .5f}, {: .5f}",
                     v7, norm(v7), angle(v1, v7), angle(v1, v7) / pi);
        fmt::println("v8 = {: .5f}, norm(v8) = {: .5f}, angle(v1,v8) = {: .5f}, {: .5f}",
                     v8, norm(v8), angle(v1, v8), angle(v1, v8) / pi);

        CHECK(std::abs(angle(v1, v1) - 0.0) < eps);
        CHECK(std::abs(angle(v1, v2) - pi * 0.25) < eps);
        CHECK(std::abs(angle(v1, v3) - pi * 0.5) < eps);
        CHECK(std::abs(angle(v1, v4) - pi * 0.75) < eps);
        CHECK(std::abs(angle(v1, v5) - pi) < eps);
    }

    TEST_CASE("Vec2d operations - wedge")
    {
        fmt::println("");
        fmt::println("Vec2d operations - wedge:");

        vec2d v1{1.0, 0.0};
        vec2d v2{normalized(vec2d(1.0, 1.0))};
        vec2d v3{0.0, 1.0};
        vec2d v4{normalized(vec2d(-1.0, 1.0))};
        vec2d v5{-1.0, 0.0};
        vec2d v6{normalized(vec2d(-1.0, -1.0))};
        vec2d v7{0.0, -1.0};
        vec2d v8{normalized(vec2d(1.0, -1.0))};

        using std::numbers::pi;

        fmt::println("v1 = {: .5f}, norm(v1) = {: .5f}, wedge(v1,v1) = {: .5f}, "
                     "sin(angle) = {: .5f}",
                     v1, norm(v1), wedge(v1, v1), sin(angle(v1, v1)));
        fmt::println("v2 = {: .5f}, norm(v2) = {: .5f}, wedge(v1,v2) = {: .5f}, "
                     "sin(angle) = {: .5f}",
                     v2, norm(v2), wedge(v1, v2), sin(angle(v1, v2)));
        fmt::println("v3 = {: .5f}, norm(v3) = {: .5f}, wedge(v1,v3) = {: .5f}, "
                     "sin(angle) = {: .5f}",
                     v3, norm(v3), wedge(v1, v3), sin(angle(v1, v3)));
        fmt::println("v4 = {: .5f}, norm(v4) = {: .5f}, wedge(v1,v4) = {: .5f}, "
                     "sin(angle) = {: .5f}",
                     v4, norm(v4), wedge(v1, v4), sin(angle(v1, v4)));
        fmt::println("v5 = {: .5f}, norm(v5) = {: .5f}, wedge(v1,v5) = {: .5f}, "
                     "sin(angle) = {: .5f}",
                     v5, norm(v5), wedge(v1, v5), sin(angle(v1, v5)));
        fmt::println("v6 = {: .5f}, norm(v6) = {: .5f}, wedge(v1,v6) = {: .5f}, "
                     "sin(angle) = {: .5f}",
                     v6, norm(v6), wedge(v1, v6), sin(angle(v1, v6)));
        fmt::println("v7 = {: .5f}, norm(v7) = {: .5f}, wedge(v1,v7) = {: .5f}, "
                     "sin(angle) = {: .5f}",
                     v7, norm(v7), wedge(v1, v7), sin(angle(v1, v7)));
        fmt::println("v8 = {: .5f}, norm(v8) = {: .5f}, wedge(v1,v8) = {: .5f}, "
                     "sin(angle) = {: .5f}",
                     v8, norm(v8), wedge(v1, v8), sin(angle(v1, v8)));

        CHECK(std::abs(wedge(v1, v1) - sin(angle(v1, v1))) < eps);
        CHECK(std::abs(wedge(v1, v2) - sin(angle(v1, v2))) < eps);
        CHECK(std::abs(wedge(v1, v3) - sin(angle(v1, v3))) < eps);
        CHECK(std::abs(wedge(v1, v4) - sin(angle(v1, v4))) < eps);
        CHECK(std::abs(wedge(v1, v5) - sin(angle(v1, v5))) < eps);
        CHECK(std::abs(wedge(v1, v6) + sin(angle(v1, v6))) < eps); // other orientation
        CHECK(std::abs(wedge(v1, v7) + sin(angle(v1, v7))) < eps); // other orientation
        CHECK(std::abs(wedge(v1, v8) + sin(angle(v1, v8))) < eps); // other orientation
    }

    TEST_CASE("Vec2d operations - project / reject")
    {
        fmt::println("");
        fmt::println("Vec2d operations - project / reject:");

        vec2d v1{5.0, 1.0};
        vec2d v2{2.0, 2.0};

        vec2d v3{project_onto(v1, v2)};
        vec2d v4{reject_from(v1, v2)};

        vec2d v5{v3 + v4};

        fmt::println("v1 = {: .5f}, norm(v1) = {: .5f}", v1, norm(v1));
        fmt::println("v2 = {: .5f}, norm(v2) = {: .5f}", v2, norm(v2));
        fmt::println("v3 = {: .5f}, norm(v3) = {: .5f}", v3, norm(v3));
        fmt::println("v4 = {: .5f}, norm(v4) = {: .5f}", v4, norm(v4));
        fmt::println("v5 = {: .5f}, norm(v5) = {: .5f}", v5, norm(v5));

        CHECK(v5 == v1);

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
}


// TEST_SUITE("fact(n):")
// {
//     TEST_CASE("fact(n): specific values")
//     {
//         CHECK(hd::fact(0) == 1.0);
//         CHECK(hd::fact(1) == 1.0);
//         CHECK(hd::fact(2) == 2.0);
//         CHECK(hd::fact(3) == 6.0);
//         CHECK(hd::fact(10) == 3628800.);
//     }
//     TEST_CASE("fact(n): throws on n<0")
//     {
//         CHECK_THROWS(hd::fact(-1));
//     }
// }