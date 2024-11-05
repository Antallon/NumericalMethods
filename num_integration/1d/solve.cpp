#include "solve.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;



// double f(double x) 
// {
//     // return 1;
//     // return x;
//     // return x*x;
//     // return x*x*x;
//     // return x*x*x*x*x;
//     return x*x*x*x*x*x*x*x*x;
//     // return x*x - 3*x + 2;
//     // return cos(x);
//     // return std::fabs(x);
//     // return 1.0/(1.0 + 25*x*x);
// }


/**
 * Simpson's rule for numerical integration.
 *
 * Given a function f(x) and a interval [a, b], this function will approximate
 * the integral of f(x) on [a, b].
 *
 * @param a The lower bound of the interval.
 * @param b The upper bound of the interval.
 * @param f The function to integrate.
 *
 * @return The approximate integral of f(x) on [a, b].
 */
double Integral_Simpson
(
    const double a, 
    const double b, 
    function<double(double)> f
)
{
    return (b - a) / 6 * (f(a) + 4 * f((a + b) / 2) + f(b));
}

/**
 * Gauss-Legendre quadrature for numerical integration.
 *
 * This function approximates the integral of a given function f(x) 
 * over the interval [a, b] using a Gauss-Legendre quadrature rule 
 * with three-point evaluation.
 *
 * @param a The lower bound of the interval.
 * @param b The upper bound of the interval.
 * @param f The function to integrate.
 *
 * @return The approximate integral of f(x) on [a, b].
 */
double Integral_Gauss
(   const double a, 
    const double b, 
    function<double(double)> f
)
{
    double x_zero = (a + b) / 2;
    double x_plus = x_zero + (b - a) / 2 * sqrt(3/5);
    double x_minus = x_zero - (b - a) / 2 * sqrt(3/5);
    return (b - a) / 18 * (5 * f(x_minus) + 8 * f(x_zero) + 5 * f(x_plus));
}

/**
 * Simpson's rule for numerical integration with composite rule.
 *
 * Given a function f(x) and a interval [a, b], this function will approximate
 * the integral of f(x) on [a, b] using a composite Simpson's rule.
 *
 * @param N The number of intervals.
 * @param a The lower bound of the interval.
 * @param b The upper bound of the interval.
 * @param f The function to integrate.
 *
 * @return The approximate integral of f(x) on [a, b].
 */
double Integral_Simpson_Composite
(
    const int N, 
    const double a, 
    const double b, 
    function<double(double)> f
)
{
    double step = (b - a) / N;
    double ans = 0;

    for(int i = 0; i < N; i++) {
        ans += Integral_Simpson(a + i * step, a + (i + 1) * step, f);
    } 

    return ans;
}

/**
 * Gauss-Kronrod's rule for numerical integration with composite rule.
 *
 * Given a function f(x) and a interval [a, b], this function will approximate
 * the integral of f(x) on [a, b] using a composite Gauss-Kronrod's rule.
 *
 * @param N The number of intervals.
 * @param a The lower bound of the interval.
 * @param b The upper bound of the interval.
 * @param f The function to integrate.
 *
 * @return The approximate integral of f(x) on [a, b].
 */
double Integral_Gauss_Composite
(
    const int N, 
    const double a, 
    const double b, 
    function<double(double)> f
)
{
    double step = (b - a) / N;
    double ans = 0;

    for(int i = 0; i < N; i++) {
        ans += Integral_Gauss(a + i * step, a + (i + 1) * step, f);
    } 

    return ans;
}



