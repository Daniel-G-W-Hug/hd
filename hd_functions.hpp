#ifndef HD_FUNCTIONS_H
#define HD_FUNCTIONS_H

#include <algorithm> // std::clamp (c++17)
#include <cmath>     // log, floor
#include <stdexcept> //std::invalid_argument

namespace hd
{
////////////////////////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////////////////////////

double linear_step(double lo_x, double hi_x, double x);
double smooth_step(double lo_x, double hi_x, double x);
double smoother_step(double lo_x, double hi_x, double x);

double log_gamma(double xx);
double fact(int n);
double log_fact(int n);
double bico(int n, int k);

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

//******************************************************************************
//                  0.                for x < lo_x
// linear_step(x) = x                 for lo_x <= x <= hi_x
//                  1.                for x > hi_x
//******************************************************************************
double linear_step(double lo_x, double hi_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    x = std::clamp((x - lo_x) / (hi_x - lo_x), 0.0, 1.0);
    return x;
}

//******************************************************************************
//                  0.                for x < lo_x
// smooth_step(x) = 3*x^2 - 2*x^3     for lo_x <= x <= hi_x
//                  1.                for x > hi_x
//
// origin: 3rd order polynomial with df/dx=0 at lo_x and hi_x
//******************************************************************************
double smooth_step(double lo_x, double hi_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    x = std::clamp((x - lo_x) / (hi_x - lo_x), 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

//******************************************************************************
//                    0.                          for x < lo_x
// smoother_step(x) = 6*x^5 - 15*x^4 + 10*x^3     for lo_x <= x <= hi_x
//                    1.                          for x > hi_x
// origin: 5th order polynomial with df/dx=0 and df/dx2=0 at lo_x and hi_x
//******************************************************************************
double smoother_step(double lo_x, double hi_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    x = std::clamp((x - lo_x) / (hi_x - lo_x), 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.);
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
    static const double cof[6] = {76.18009172947146, -86.50532032941677,
                                  24.01409824083091, -1.231739572450155,
                                  0.1208650973866179e-2, -0.5395239384953e-5};
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * std::log(tmp);
    ser = 1.000000000190015;

    for (int j = 0; j <= 5; ++j)
    {
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

    if (n < 0)
    {
        throw std::invalid_argument("Negative argument in fact(n).");
    }

    if (n > 32)
    {
        // might overflow depending on n and sizeof(double)
        return std::exp(log_gamma(n + 1.0));
    }

    while (ntop < n)
    {
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

    if (n < 0)
    {
        throw std::invalid_argument("Negative argument in log_fact(n).");
    }

    if (n <= 1)
    {
        return 0.0;
    }

    if (n <= 100)
    {
        // in range of table
        return a[n] ? a[n] : (a[n] = log_gamma(n + 1.0));
    }
    else
    {
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

} // namespace hd

#endif // HD_FUNCTIONS_H