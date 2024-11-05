#include "solve.hpp"
#include "polynom.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;

/**
 * @brief Computes the value of Lagrange polynomial at a given point
 * 
 * Computes the value of Lagrange polynomial at a given point x using n nodes
 * and corresponding function values y.
 * 
 * @param[in] n Number of nodes
 * @param[in] x Point at which to compute the polynomial
 * @param[in] nodes Vector of length n containing the nodes
 * @param[in] y Vector of length n containing the function values
 * @return The value of the Lagrange polynomial at x
 */
double L
(
    const int n, 
    const double x, 
    const vector<double>& nodes, 
    const vector<double>& y
)
{

    double sum = 0.0;
    for(int i=0; i < n; i++)
    {
        double mult = 1.0;
        for(int j=0; j<n; j++)
        {
            if(j==i) continue;
            else
            {
                mult *= (x - nodes[j])/(nodes[i] - nodes[j]);
            }
            

        }
        sum += mult * y[i];

    }
    return sum;
}

/**
 * @brief Computes the value of a polynomial at a given point
 * 
 * Computes the value of a polynomial of degree n-1 at a given point x using
 * the coefficients in the vector coeff.
 * 
 * @param[in] n Number of coefficients.
 * @param[in] x Point at which to compute the polynomial.
 * @param[in] coeff Vector containing the coefficients of the polynomial.
 * 
 * @return The value of the polynomial at x.
 */
double P
(
    const int n, 
    const double x, 
    const vector<double>& coeff
)
{
    double sum = 0.0;
    for(int i = 0; i < n; i++)
    {  
        sum += coeff[i] * pow(x,i);
    }
    return sum;    

}



/**
 * @brief Fills a matrix with powers of nodes.
 * 
 * This function populates the matrix A with powers of each element in the
 * nodes vector. Each row i of the matrix corresponds to the powers of 
 * nodes[i], starting from 0 up to n-1.
 * 
 * @param[in] n Number of nodes and the degree of the polynomial.
 * @param[out] A Vector representing the matrix to be filled. It must have
 * a size of n*n.
 * @param[in] nodes Vector containing the nodes for which the matrix is
 * populated with powers.
 */
void Pmatfull
(
    const int n, 
    vector<double>& A, 
    const vector<double>& nodes
)
{
    for(int i=0; i<n; i++)
    {
        double p = 1.0;
        for(int j=0; j<n; j++)
        {
            
            A[i*n + j] = p;
            p *= nodes[i];
        }
    }
}