#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <vector>

// include functions to be tests
#include "hd_ga.hpp"

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

    TEST_CASE("Vec2d operations")
    {
        fmt::println("");
        fmt::println("Vec2d operations:");

        vec2d v1{2.0, 1.0};
        vec2d v2{normalized(v1)};

        vec2d v3{2.0, 6.0};
        vec2d v4{inverse(v3)};

        fmt::println("v1 = {}, norm(v1) = {}", v1, norm(v1));
        fmt::println("v2 = normalized(v1) = {}, norm(v2) = {}", v2, norm(v2));

        CHECK(std::abs(sq_norm(v1) - 5.0) < eps); // comparison (equality)
        CHECK(std::abs(sq_norm(v2) - 1.0) < eps); // comparison (inequality)
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);

        fmt::println("v3 = {}, norm(v3) = {}", v3, norm(v3));
        fmt::println("v4 = inverse(v3) = {}, norm(v4) = {}", v4, norm(v4));

        vec2d v5{3.0, 0.0};
        vec2d v6{2.0, 2.0};
        fmt::println("v5 = {}, norm(v5) = {}", v5, norm(v5));
        fmt::println("v6 = {}, norm(v6) = {}", v6, norm(v6));
        fmt::println("dot(v5,v6) = {}", dot(v5, v6));
        fmt::println("project_onto(v6,v5) = {}", project_onto(v6, v5));
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