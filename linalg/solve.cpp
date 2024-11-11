#include "solve.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;




/**
 * @brief Fill matrix A with coefficients for finite difference method
 * 
 * @param[in] N - size of matrix A
 * @param[in] p - coefficient from equation
 * @param[out] A - matrix A with coefficients for finite difference method
 */
void fillmat
(
    const int N, 
    const double p, 
    vector<double>& A 
)
{
    for(int i = 0; i < N + 1; i++)
    {
        for(int j = 0; j < N + 1; j++)
        {
            if(i == j && i > 0 && i < N && j > 0 && j < N)
            A[i * (N + 1) + j] = p + 2.0 * (double)(N) * (double)(N);

            else if((std::abs(i - j) == 1) && i > 0 && i < N && j > 0 && j < N) 
            A[i * (N + 1) + j] = -1.0 * (double)(N) * (double)(N);
            
            else A[i * (N + 1) + j] = 0.0;
        }
    }
}



/**
 * @brief Fill vector b with the right side of the system of equations
 * 
 * @param[in] A - matrix A with coefficients for finite difference method
 * @param[in] N - size of matrix A
 * @param[out] b - vector b with the right side of the system of equations
 */
void fillf
(
    const vector<double>& A, 
    const int N, 
    vector<double>& f
)
{
    f[0] = 0.0;
    f[N] = 0.0;
    for (int k = 1; k < N; k++)
    {
        f[k] = 0;
        for (int j = 0; j < N; j += 2)
        {
            f[k] += A[k * N + j];
        }
    }
}


/**
 * @brief Multiply matrix A by vector b
 * 
 * @param[in] N - size of matrix A
 * @param[in] A - matrix A with coefficients for finite difference method
 * @param[in] b - vector b with the right side of the system of equations
 * @param[out] res - vector res with result of multiplication
 */
void multiply
(
    const int N,
    const vector<double>& A, 
    const vector<double>& b, 
    vector<double>& res
) 
{
    for(int i = 0; i < N + 1; i++)
    {
        res[i] = 0.0;
        for(int j = 0; j < N + 1; j++)
        {
            res[i] += A[i * (N + 1) + j] * b[j];
        }
    } 
}



double vecnorm(const vector<double>& v) 
{
    double sum = 0.0;
    for(size_t i = 0; i < v.size(); i++)
    {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

void vecdiff(const vector<double>& v1, const vector<double>& v2, vector<double>& diff)
{
    for(size_t i = 0; i < diff.size(); i++)
    {
        diff[i] = v1[i] - v2[i];
    }
}

/**
 * @brief Compute eigen vector of the matrix A
 * 
 * @param[in] N - size of matrix A
 * @param[in] m - number of eigen vector
 * @param[in] k - number of element in eigen vector
 * @return eigen vector of the matrix A
 */
double compute_eigen_vector
(
    const int N, 
    const int m, 
    const int k 
)
{
    return sin(M_PI * m * k / (double)(N));
}



/**
 * @brief Compute eigen value of the matrix A
 * 
 * @param[in] N - size of matrix A
 * @param[in] m - number of eigen value
 * @param[in] p - parameter of the problem
 * @return eigen value of the matrix A
 */
double compute_eigen_value
(
    const int N, 
    const int m, 
    const double p
)
{
    double l = p - 2 * N * N * (cos(M_PI * m / (double)(N)) - 1);
    return l;
}




/**
 * @brief Compute the coefficient d for the finite difference method
 * 
 * @param[in] N - size of the problem
 * @param[in] m - index of the eigen vector
 * @param[in] p - parameter of the problem (unused in this function)
 * @param[in] f - vector of function values
 * @return coefficient d based on the given eigen vector and function values
 */
double compute_coeff_d
(
    const int N, 
    const int m, 
    const double p, 
    const vector<double>& f
)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        sum += 2 * f[i] * compute_eigen_vector(N, m, i) / (double)(N);
       
    }

    return sum;
}




/**
 * @brief Compute the solution vector y for the problem
 * 
 * @param[in] N - size of the problem
 * @param[in] p - parameter of the problem
 * @param[out] y - solution vector
 * @param[in] f - vector of function values
 */
void compute_solution_y
(
    const int N, 
    const double p,
    const vector<double>& f,
    vector<double>& y

)
{
    for(int k = 0; k < N + 1; k++)
    {
        y[k] = 0.0;
        for(int m = 1; m < N; m++)
        {
           
            y[k] += compute_coeff_d(N, m, p, f) /compute_eigen_value(N, m, p) * compute_eigen_vector(N, m, k) ;
        }
    }

}


void Richardson
(
    vector<double>&x, 
    const vector<double>& A, 
    const vector<double>& b, 
    const double tau, 
    const int N, 
    const int mIter,
    vector<double>& prod
)
{
    double norma = 0.0;

    for(int i = 0; i < N + 1; i++)
    {
        x[i]    = 0.0;
        prod[i] = 0.0;
    }

    for(int m = 0; m < mIter; m++)
    {
        multiply(N, A, x, prod);
        for(int i = 0; i < N + 1; i++)
        {
            x[i] = x[i] - tau * prod[i] + tau * b[i];
        }
    }

    // cout<<"Iter number = "<<mIter<<endl;
    // for(int i = 0; i < N + 1; i++) cout<< x[i];
    // cout<<endl<<endl;


    // for(int i = 0; i < N; i++) prod[i] = 0.0;

    // multiply(N, A, x, prod);
    // for(int i = 0; i < N; i++)
    // {
    //     norma += (prod[i] - b[i]) * (prod[i] - b[i]);
    // }
    // return sqrt(norma);

}











void printMatrix(vector<double>& matrix,int n, int m,  vector<double>& vector_b){
        int l;
        if (n <= 10)
                l = n;
        else
                l = 5;
        for (int i = 0; i < m; i++){
                for (int j = 0; j < l; j++){
                        if ((i*n + j) <= (n*m - 1))
                                printf("%10.3e ", matrix[i*n + j]);
                }
                printf("  |    %10.3e\n", vector_b[i]);
        }
        printf("\n");
}