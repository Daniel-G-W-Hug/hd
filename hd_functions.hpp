#ifndef HD_FUNCTIONS_H
#define HD_FUNCTIONS_H

#include <algorithm> // std::clamp (c++17)

namespace hd
{

double linear_step(double lo_x, double hi_x, double x);
double smooth_step(double lo_x, double hi_x, double x);
double smoother_step(double lo_x, double hi_x, double x);

/*****************************************************************************
//                  0.                for x < lo_x
// linear_step(x) = x                 for lo_x <= x <= hi_x
//                  1.                for x > hi_x
******************************************************************************/
double linear_step(double lo_x, double hi_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    x = std::clamp((x - lo_x) / (hi_x - lo_x), 0.0, 1.0);
    return x;
}

/*****************************************************************************
//                  0.                for x < lo_x
// smooth_step(x) = 3*x^2 - 2*x^3     for lo_x <= x <= hi_x
//                  1.                for x > hi_x
//
// origin: 3rd order polynomial with df/dx=0 at lo_x and hi_x
******************************************************************************/
double smooth_step(double lo_x, double hi_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    x = std::clamp((x - lo_x) / (hi_x - lo_x), 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

/*****************************************************************************
//                    0.                          for x < lo_x
// smoother_step(x) = 6*x^5 - 15*x^4 + 10*x^3     for lo_x <= x <= hi_x
//                    1.                          for x > hi_x
// origin: 5th order polynomial with df/dx=0 and df/dx2=0 at lo_x and hi_x
******************************************************************************/
double smoother_step(double lo_x, double hi_x, double x)
{
    // scale, bias and saturate x to 0..1 range
    x = std::clamp((x - lo_x) / (hi_x - lo_x), 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.);
}

} // namespace hd

#endif // HD_FUNCTIONS_H