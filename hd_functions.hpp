#ifndef HD_FUNCTIONS_H
#define HD_FUNCTIONS_H

#include <algorithm> // std::clamp (c++17)
#include <cmath>     // log, floor
#include <limits>    // numerical limits
#include <stdexcept> //std::invalid_argument

namespace hd {
////////////////////////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////////////////////////

double linear_step(double low_x, double high_x, double x);
double smooth_step(double low_x, double high_x, double x);
double smoother_step(double low_x, double high_x, double x);

double log_gamma(double xx);
double fact(int n);
double log_fact(int n);
double bico(int n, int k);

enum class split_t { geometric, arithmetic };
int oo_magnitude(double x, split_t s = split_t::geometric);

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

//******************************************************************************
//                  0.                               for x < low_x
// linear_step(x) = (x - low_x)/(high_x - low_x)     for low_x <= x <= high_x
//                  1.                               for x > high_x
//******************************************************************************
double linear_step(double low_x, double high_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    x = std::clamp((x - low_x) / (high_x - low_x), 0.0, 1.0);
    return x;
}

//******************************************************************************
//                  0.                for x < low_x
// smooth_step(x) = 3*x^2 - 2*x^3     for low_x <= x <= high_x
//                  1.                for x > high_x
//
// origin: 3rd order polynomial with df/dx=0 at low_x and high_x
//******************************************************************************
double smooth_step(double low_x, double high_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    // with  x = (x-low_x)/(high_x-low_x)    normalize x to range (0.,1.)
    x = std::clamp((x - low_x) / (high_x - low_x), 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

//******************************************************************************
//                    0.                          for x < low_x
//          with  x = (x-low_x)/(high_x-low_x)    normalize x to range (0.,1.)
// smoother_step(x) = 6*x^5 - 15*x^4 + 10*x^3     for low_x <= x <= high_x
//                    1.                          for x > high_x
// origin: 5th order polynomial with df/dx=0 and df/dx2=0 at low_x and high_x
//******************************************************************************
double smoother_step(double low_x, double high_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    x = std::clamp((x - low_x) / (high_x - low_x), 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

////////////////////////////////////////////////////////////////////////////////
// functions as per "Numerical recipes in C", 2nd edition, Chapter 6
//
// n! = gamma(n+1)
//
// with gamma() expressed as log_gamma() to avoid overflow
////////////////////////////////////////////////////////////////////////////////

//******************************************************************************
// returns log(gamma(xx)) für xx>0.
//******************************************************************************
double log_gamma(double xx)
{
    double x, y, tmp, ser;
    static const double cof[6] = {76.18009172947146,     -86.50532032941677,
                                  24.01409824083091,     -1.231739572450155,
                                  0.1208650973866179e-2, -0.5395239384953e-5};
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * std::log(tmp);
    ser = 1.000000000190015;

    for (size_t j = 0; j <= 5; ++j) {
        ser += cof[j] / ++y;
    }
    return -tmp + std::log(2.5066282746310005 * ser / x);
}

//******************************************************************************
// return n! as double value
//******************************************************************************
double fact(int n)
{
    static int ntop = 6;

    // fill in only as required (static arrays are zero initialized)
    static double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0};

    if (n < 0) {
        throw std::invalid_argument("Negative argument in fact(n).");
    }

    if (n > 32) {
        // might overflow depending on n and sizeof(double)
        return std::exp(log_gamma(n + 1.0));
    }

    while (ntop < n) {
        // fill our array further if needed, for potential reuse in static array
        int j = ntop++;
        a[ntop] = a[j] * ntop;
    }
    return a[n];
}

//******************************************************************************
// returns log(n!)
//******************************************************************************
double log_fact(int n)
{
    // static arrays are zero initialized
    static double a[101];

    if (n < 0) {
        throw std::invalid_argument("Negative argument in log_fact(n).");
    }

    if (n <= 1) {
        return 0.0;
    }

    if (n <= 100) {
        // in range of table
        return a[n] ? a[n] : (a[n] = log_gamma(n + 1.0));
    }
    else {
        // out of range of table
        return log_gamma(n + 1.0);
    }
}

//******************************************************************************
// returns binomial coefficient (n choose k) = n!/( k!*(n-k)! ) for 0 <= k <= n
//******************************************************************************
double bico(int n, int k)
{
    // std::floor is for clean up of roundoff errors for smaller values of n and k
    return std::floor(0.5 + std::exp(log_fact(n) - log_fact(k) - log_fact(n - k)));
}

//******************************************************************************
// returns the order of magnitude (oom) of a given number as exponent 10^N
//
//
// geometric split: with N such that 10^(-1/2) <= abs(x)/10^N < 10^(1/2)
//
//   => 10^(N-1/2) <= abs(x) < 10^(N+1/2) with 10^N as geometric midpoint
//
//
// arithmetic split: with N such that 0.5 < abs(x)/10^N <= 5.0
//
//   => 0.5*10^N < abs(x) <= 5.0*10^N
//
//******************************************************************************

int oo_magnitude(double x, split_t s)
{

    if (abs(x) <= std::numeric_limits<double>::min()) return 0;

    if (s == split_t::geometric) // geometric split (default)
        return static_cast<int>(std::floor(std::log10(std::abs(x)) + 0.5));
    else if (s == split_t::arithmetic) //  arithmetric split
        return static_cast<int>(std::log10(std::abs(x) / 0.5));
    else throw std::invalid_argument("hd::oo_magnitude: case in split_t not handled.");
}

// Kronecker delta function
template <typename T = double>
    requires std::same_as<T, double> || std::same_as<T, float> || std::same_as<T, int> ||
             std::same_as<T, long> || std::same_as<T, long long> ||
             std::same_as<T, std::size_t>
constexpr T kronecker(size_t i, size_t j)
{
    return (i == j) ? T(1) : T(0);
}

// concept to restrict templates to int or size_t
template <typename T>
concept IntegerIndex = std::same_as<T, int> || std::same_as<T, std::size_t>;

// eps is the Levi-Civita symbol (permutation symbol) for n indices
// eps(i,j,...) = +1 for even permutations
//              = -1 for odd permutations
//              =  0 if any index is repeated
template <IntegerIndex... Args> auto eps(Args... args)
{
    constexpr std::size_t n = sizeof...(Args);

    // store arguments in arrary
    int indices[] = {static_cast<int>(args)...};

#ifndef _HD_SKIP_EPS_INDEX_RANGE_TEST
    // check if indices are a permissible permutation
    int sorted[n];
    std::copy(std::begin(indices), std::end(indices), sorted);
    std::sort(std::begin(sorted), std::end(sorted));

    // Either {0,1,2,...,n-1} or {1,2,3,...,n} is permissible

    // check for duplicate indices first
    bool has_duplicates = false;
    for (std::size_t i = 0; i < n - 1; ++i) {
        if (sorted[i] == sorted[i + 1]) {
            has_duplicates = true;
            break;
        }
    }

    // only if there are NO duplicates, check if permissible permutation
    if (!has_duplicates) {
        bool valid_from_zero = true;
        bool valid_from_one = true;

        for (std::size_t i = 0; i < n; ++i) {
            if (sorted[i] != static_cast<int>(i)) valid_from_zero = false;
            if (sorted[i] != static_cast<int>(i + 1)) valid_from_one = false;
        }

        if (!valid_from_zero && !valid_from_one) {
            throw std::invalid_argument("eps: Indices must be a permutation of "
                                        "consecutive integers starting from 0 or 1");
        }
    }
#endif

    // product of all pair differences for (j > i)
    int numerator = 1;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            numerator *= (indices[j] - indices[i]);
        }
    }

    // product of all differences of positions (j > i)
    int denominator = 1;
    for (int i = 0; i < static_cast<int>(n); ++i) {
        for (int j = i + 1; j < static_cast<int>(n); ++j) {
            denominator *= (j - i);
        }
    }

    // normalize to +1 or -1
    return numerator / denominator;
}

} // namespace hd

#endif // HD_FUNCTIONS_H