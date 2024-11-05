#include <math.h>
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
using namespace std;



double u(double x, double y)
{
    // return cos(2 * M_PI * x)*cos(2 * M_PI * y);
    // return sin(2 * M_PI * x)*sin(2 * M_PI * y);
    return cos(2 * M_PI * x)*sin(2 * M_PI * y);
    // return cos(2 * M_PI * x)*cos(2 * M_PI * y) * cos(2 * M_PI * x)*cos(2 * M_PI * y);


    // return x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y);
    // return x*x*(1-x*x) * sin(M_PI * x)*y*y*(1-y*y) * sin(M_PI * y);
    // return x*y;
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
 * @brief Compute trigonometric functions for given x and N
 * @param[in] N Number of trigonometric functions
 * @param[in] x Value for which to compute trigonometric functions
 * @param[out] cx Vector of computed trigonometric functions
 */
void computeTrig(int N, double x, vector<double>& cx)
{
    for(int m=0; m < N; m++)
    {
        cx[m] = cos(M_PI * m * x);
    }
}



/**
 * @brief Compute zero order moment of a given 2D array
 * @param[in] N Number of points in the array
 * @param[in] u 2D array of values
 * @return Zero order moment of the array
 */
double compute_zero(int N, vector<double>& u)
{
    double h = 1.0 / ( N - 0.5 );
    double zero = 0.0;
    for (int i = 0; i < N; i++)
    {
        for(int j=0; j < N; j++)
        {
            zero += u[i * N +j] * h * h ;
        }
    }

    return zero; 
}


/**
 * @brief Compute coefficients for a given function using trigonometric basis functions.
 * 
 * @param[in] N Number of points for computation.
 * @param[in] x Vector of grid points.
 * @param[out] coef Iterator to the beginning of the output coefficients.
 * @param[in] u Iterator to the input function values.
 * 
 * This function calculates the coefficients of a given function
 * represented by its values at specified grid points using a trigonometric
 * basis. The coefficients are computed by integrating the product of the
 * function values and the cosine of the basis functions over the grid points.
 */
void coeff(int N , vector<double>& x, vector<double>::iterator coef, vector<double>::iterator u)
{
    double h = 1.0 / ( N - 0.5 );
    
    for(int i=0; i < N; i++)
    {
        *(coef + i) = 0.0;
        for(int j=1; j < N; j++)
        {
            *(coef + i) += *(u+j) * h * cos(M_PI * i * x[j]);
        }
        if(i!=0) *(coef + i) *= 2;
    }

}

/**
 * @brief Fill a 2D array with values of a given function at specified grid points
 * 
 * @param[in] N Number of points in the array
 * @param[out] U 2D array to fill with values of the function
 * @param[in] points Vector of grid points
 * 
 * This function computes the values of a given function
 * at specified grid points and stores them in a 2D array.
 */
void Ufull(int N, vector<double>& U, vector<double>& points)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            
            U[i * N + j] =  u(points[i], points[j]);
            
        }
    }
}


/**
 * @brief Compute coefficients for a given 2D array using trigonometric basis functions.
 * 
 * This function calculates the coefficients of a given 2D array represented by
 * its values at specified grid points using a trigonometric basis. The coefficients
 * are computed by integrating the product of the function values and the cosine of
 * the basis functions over the grid points.
 * 
 * @param[in] N Number of points in the array
 * @param[in] x Vector of grid points
 * @param[out] C 2D array to store the computed coefficients
 * @param[in] U 2D array of values to compute coefficients for
 */
void Cfull(int N, vector<double>& x,  vector<double>& C, vector<double>& U)
{
  

    for(int i = 0; i < N; i++)
    {
        // C[i*N] = 0;
        coeff(N,x,C.begin() + i * N, U.begin() + i * N);
    }
}

/**
 * @brief Compute coefficients for a given 2D array transposed using trigonometric basis functions.
 * 
 * This function calculates the coefficients of a given 2D array represented by
 * its values at specified grid points using a trigonometric basis. The coefficients
 * are computed by integrating the product of the function values and the cosine of
 * the basis functions over the grid points. The input array is transposed, so that
 * instead of computing the coefficients for each row, the function computes the
 * coefficients for each column of the input array.
 * 
 * @param[in] N Number of points in the array
 * @param[in] x Vector of grid points
 * @param[in] C 2D array of values to compute coefficients for
 * @param[out] D 2D array to store the computed coefficients
 */
void Dfull(int N, vector<double>& x,  vector<double>& C, vector<double>& D)
{
    for(int i = 0; i < N; i++)
    {
        vector <double> col(N);
        // D[i*N] = 0;
        for(int h = 0; h < N; h++)
        {
            col[h] = C[h * N + i];
        }
        coeff(N, x, D.begin() + i * N, col.begin());
    }
}


/**
 * @brief Computes an approximation of a 2D function at a given point (x, y) using Fourier coefficients.
 * 
 * This function calculates an approximation of a 2D function's value at a specified 
 * point (x, y) using precomputed Fourier coefficients stored in the vector `D`. It 
 * utilizes trigonometric basis functions computed by the `computeTrig` function for 
 * both x and y coordinates. The approximation is obtained by summing the product of 
 * the coefficients and the trigonometric values.
 * 
 * @param[in] N Number of terms in the Fourier series.
 * @param[in] x The x-coordinate of the point at which to compute the approximation.
 * @param[in] y The y-coordinate of the point at which to compute the approximation.
 * @param[in] points Vector containing the grid points.
 * @param[in] U Vector containing values used for computing the approximation (not used in this function).
 * @param[in] D Vector containing the precomputed Fourier coefficients.
 * 
 * @return The computed approximation of the 2D function at the point (x, y).
 */
double aprox2 (int N, double x, double y, vector<double>& points, vector<double>& U, vector<double>& D)
{
    vector <double> cx(N);
    vector <double> cy(N);
    computeTrig(N, x, cx);
    computeTrig(N, y, cy);
    // double z = compute_zero(N, U);
    // cout<<z<<endl;
    double z = 0.0;
    for(int i=0; i < N; i++)
    {
        for(int j=0; j < N; j++)
        {
            z += D[i * N + j] * cx[i] * cy[j];
        }
    }
    return z;
}








