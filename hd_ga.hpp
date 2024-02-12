#ifndef HD_GA_HPP
#define HD_GA_HPP

// author: Daniel Hug, 2024

////////////////////////////////////////////////////////////////////////////////
// provide data types for representing GA in 2D and 3D
// e.g. scalar, vector, bivector, trivector
// and corresponding operations
////////////////////////////////////////////////////////////////////////////////

#include <array>
#include <cmath>    // sqrt, abs
#include <compare>  // <=>
#include <concepts> // std::floating_point<T>
#include <cstdint>  // uint8_t
#include <iostream>

namespace hd {

////////////////////////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////////////////////////

// create an algebra G(P,N,Z) with:
// P generators for numbers that square to +1
// N generators for numbers that square to -1
// Z generators for numbers that square to  0
//
// used to provide frequently used values of the algebra
// should be assigned to a const variable
template <uint8_t P, uint8_t N = 0, uint8_t Z = 0>
    requires(P + N + Z >= 2 && P + N + Z <= 4 && P >= 2 && P <= 4 && N == 0 &&
             Z == 0) // no implementation for other algebras yet
struct Algebra {
    constexpr static const uint8_t p() { return P; }; // number of +1 generators
    constexpr static const uint8_t n() { return N; }; // number of -1 generators
    constexpr static const uint8_t z() { return Z; }; // number of  0 generators
    constexpr static const uint8_t dim_space()
    {
        return P + N + Z;
    }; // dimension of the space
    constexpr static const uint8_t num_components()
    {
        return 1 << (P + N + Z); // the number of basis components == 2^dim_space
    };

    constexpr static const std::array<const uint8_t, dim_space() + 1>
        num_components_grade = []() -> std::array<const uint8_t, dim_space() + 1> {
        if constexpr (dim_space() == 2) return {1, 2, 1};
        if constexpr (dim_space() == 3) return {1, 3, 3, 1};
        if constexpr (dim_space() == 4) return {1, 4, 6, 4, 1};
    }();

    constexpr static const std::array<const char[6], num_components()> basis_name =
        []() -> std::array<const char[6], num_components()> {
        if constexpr (dim_space() == 2) return {"    1", "   e1", "   e2", "  e12"};
        if constexpr (dim_space() == 3)
            return {"    1", "   e1", "   e2", "   e3",
                    "  e23", "  e31", "  e12", " e123"};
        if constexpr (dim_space() == 4)
            return {"    1", "   e1", "   e2", "   e3", "   e4", "  e41",
                    "  e42", "  e43", "  e23", "  e31", "  e12", " e423",
                    " e431", " e412", " e123", "e0123"};
    }();
};


////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

template <typename T = float>
    requires(std::floating_point<T>)
struct Vec2d {

    // ctor
    Vec2d<T>() = default;
    Vec2d<T>(T x_in, T y_in) : x(x_in), y(y_in) {}

    T x{};
    T y{};

    // equality
    bool operator==(const Vec2d<T>& rhs) const
    {
        // comparison using SquaredNorm (metric not considered!)
        auto snl = x * x + y * y;
        auto snr = rhs.x * rhs.x + rhs.y * rhs.y;
        if (std::abs(snl - snr) < 1.e-7) return true;
        return false;
    }

    // comparison via operator <=> using SquaredNorm of vector for comparison
    std::strong_ordering operator<=>(const Vec2d<T>& rhs) const
    {
        // comparison using SquaredNorm (metric not considered!)
        auto snl = x * x + y * y;
        auto snr = rhs.x * rhs.x + rhs.y * rhs.y;
        if (std::abs(snl - snr) < 1.e-7) return std::strong_ordering::equal;
        if (snl < snr) return std::strong_ordering::less;
        return std::strong_ordering::greater;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Vec2d<U>& v);
};

// for printing via iostream
template <typename T> std::ostream& operator<<(std::ostream& os, const Vec2d<T>& v)
{
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}


using value_t = float;
using scalar_t = value_t;
using pseudoscalar_t = value_t;
using Vector2d = Vec2d<value_t>;

// unary minus
inline Vector2d operator-(const Vector2d& v) { return Vector2d(-v.x, -v.y); }

// adding vectors
inline Vector2d operator+(const Vector2d& v1, const Vector2d& v2)
{
    return Vector2d(v1.x + v2.x, v1.y + v2.y);
}

// substracting vectors
inline Vector2d operator-(const Vector2d& v1, const Vector2d& v2)
{
    return Vector2d(v1.x - v2.x, v1.y - v2.y);
}

// multiply a vector with a scalar
inline Vector2d operator*(const Vector2d& v, value_t s)
{
    return Vector2d(v.x * s, v.y * s);
}
inline Vector2d operator*(value_t s, const Vector2d& v)
{
    return Vector2d(v.x * s, v.y * s);
}

// devide a vector by a scalar
inline Vector2d operator/(const Vector2d& v, value_t s)
{
    value_t inv = 1.0 / s; // for multiplicaton with inverse value
    return Vector2d(v.x * inv, v.y * inv);
}


// return magnitude of vector
// (TODO: extend for non-orthonormal systems using a metric)
inline value_t Norm(const Vector2d& v) { return std::sqrt(v.x * v.x + v.y * v.y); }

// return squared magnitude of vector
// (TODO: extend for non-orthonormal systems using a metric)
inline value_t SquaredNorm(const Vector2d& v) { return v.x * v.x + v.y * v.y; }

// return a vector nomalized to Norm(v) == 1.0
inline Vector2d Normalize(const Vector2d& v)
{
    value_t inv = 1.0 / Norm(v); // for multiplication with inverse of Norm
    return Vector2d(v.x * inv, v.y * inv);
}

// return dot-product of two vectors
// (TODO: extend for non-orthonormal systems by using metric)
inline value_t Dot(const Vector2d& v1, const Vector2d& v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

// wedge product (returns a bivector, which is the pseudoscalar in 2d)
// (TODO: extend for non-orthonormal systems by using metric)
inline pseudoscalar_t Wedge(const Vector2d& v1, const Vector2d& v2)
{
    return v1.x * v2.y - v1.y * v2.x;
}

// projection of v1 onto v2 (v2 must be normalized to Norm(v2) == 1)
inline Vector2d project_onto(const Vector2d& v1, const Vector2d& v2)
{
    return v2 * Dot(v1, v2);
}

// rejection of v1 from v2 (v2 must be normalized to Norm(v2) == 1)
inline Vector2d reject_from(const Vector2d& v1, const Vector2d& v2)
{
    return v1 - v2 * Dot(v1, v2);
}


struct Frame {
    Frame(const Vector2d& v1, const Vector2d& v2) : e1(v1), e2(v2) {}
    Vector2d e1;
    Vector2d e2;
};

// T& s{&c[0]};   // scalar part
// T& x{&c[1]};   // vector part, basis vector e1
// T& y{&c[2]};   // vector part, basis vector e2
// T& pss{&c[3]}; // pseudo scalar part

namespace ega2d { // euclidean geometric algebra 2d

} // namespace ega2d

} // namespace hd

#endif // HD_GA_HPP