#include <math.h>
#include <cmath>
#include <functional>
#include <vector>
using namespace std;

double u(double x)
{
    // return cos(2* M_PI * x)/2 + cos(M_PI * x);
    // return x*x*(1-x*x) * sin(M_PI * x); 
    // return cos(M_PI * x) * cos(M_PI * x); 
    return cos(2*M_PI * x) * sin(2*M_PI * x); 

    // return x*x*(1-x)*(1-x) ;
    // + cos(M_PI * x);
    // return cos(2 * M_PI * x);
}


/**
 * @brief Compute grid points and points for numerical integration
 * @param[in] N Number of points
 * @param[out] vec Vector of grid points
 * @param[out] p Vector of points for numerical integration
 */
void computePoints(int N, vector<double>& vec,vector<double>& p)
{
    double h = 1.0 / ( N - 0.5 );
    double x_k = -h / 2;
    vec[0] = x_k;
    p[0] = x_k;

    for(int i=1; i < N; i++)
    {
        x_k += h;
        vec[i] = x_k;
    }

    for(int j = 1; j < (N-1)*3 + 1; j++)
    {
        p[j] = p[j-1] + h/3;
    }
}


/**
 * @brief Computes the trigonometric cosine values for a given x and stores them in a vector.
 * 
 * This function calculates the cosine values for multiples of a given angle (M_PI * m * x) 
 * and stores them in the provided vector `sx`. The calculation is performed for `m` ranging 
 * from 0 to N-1.
 * 
 * @param[in] N The number of trigonometric values to compute.
 * @param[in] x The angle multiplier used in the cosine computation.
 * @param[out] sx Vector to store the computed cosine values.
 */
void computeTrig(int N, double x, vector<double>& sx)
{
    // sx[0] = 0;
    for(int m=0; m < N; m++)
    {
        // sx[m] = sin(M_PI * m * x);
        sx[m] = cos(M_PI * m * x);
        // sx[m] = cos(2 * M_PI * m * x / (2*N-1));
    }
}


/**
 * @brief Simple rectangular integration algorithm.
 * 
 * This function performs a simple numerical integration by approximating the
 * integral of a given function `f` on a set of points `x` by summing the
 * function values at each point multiplied by the step size `h`.
 * 
 * @param[in] N The number of points in the `x` array.
 * @param[in] x The array of points at which to evaluate the function.
 * @param[in] f The function to integrate.
 * 
 * @return The approximate integral value.
 */
double integrate(int N, vector<double>& x, function<double(double)> f)
{
    double h = 1.0 / ( N - 0.5 );
    double integral = 0.0;
    for(int i=1; i < N; i++)
    {
        integral+= f(x[i]) * h;
    }
    return integral;
}

/**
 * @brief Computes the integral of the function u over a set of points.
 * 
 * This function calculates the integral of the function `u` across a 
 * specified set of points `x` using a simple integration technique. 
 * It utilizes the `integrate` function to perform the numerical 
 * integration and returns the computed result.
 * 
 * @param[in] N Number of points in the vector `x`.
 * @param[in] x Vector containing the points over which to integrate.
 * 
 * @return The computed integral value, representing a form of the 
 * "zero" of the function `u` over the points.
 */
double compute_zero(int N, vector<double>& x)
{
    double zero =  integrate(N ,x , [](double a){
        return u(a);
    });
    return zero;
}


/**
 * @brief Computes the Fourier coefficients for the function u.
 * 
 * This function calculates the Fourier coefficients for the function `u`
 * using a numerical integration method. It iterates through the range of
 * N-1 and computes each coefficient by integrating the product of `u` and
 * the cosine function over the specified points `x`. The computed coefficients
 * are stored in the vector `coef`.
 * 
 * @param[in] N Number of terms in the Fourier series.
 * @param[in] x Vector containing the points over which to integrate.
 * @param[out] coef Vector to store the computed Fourier coefficients.
 */
void coeff(int N , vector<double>& x, vector<double>& coef)
{

    for(int i=1; i < N; i++)
    {
        coef[i] = 2 * integrate(N, x, [i](double a){
            return u(a) * cos(M_PI * i * a);
        });
    }

}

/**
 * @brief Computes an approximation of the function u at the point x_k using
 *        the Fourier coefficients computed by the function coeff.
 *
 * @param[in] N   Number of terms in the Fourier series.
 * @param[in] x_k The point at which to compute the approximation.
 * @param[in] x   Vector containing the points used to compute the
 *                Fourier coefficients.
 * @param[in] coef Vector containing the computed Fourier coefficients.
 *
 * @return The computed approximation of the function u at the point x_k.
 */
double aprox(int N, double x_k, vector<double>& x, vector<double>& coef)
{
    vector<double> sx(N);
    computeTrig(N, x_k, sx);
    double y_k = compute_zero(N,x);
    // double y_k = 0.0;
    for(int i=1; i < N; i++)
    {
        y_k+=coef[i] * sx[i];
    }
    return y_k;
}








