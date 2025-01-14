#include "solve2.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833673362440656643086021394946395224737190702179860943 

double u(double x)
{
    // return cos(2* M_PI * x)/2 + cos(M_PI * x);
    // return x*x*(1-x*x) * sin(M_PI * x); 
    // return cos(M_PI * x) * cos(M_PI * x); 
    // return cos(2*M_PI * x) * sin(2*M_PI * x); 
    // return x*x*(1-x)*(1-x) ;
    // + cos(M_PI * x);
    // return cos(2 * M_PI * x);
    return cos(2*M_PI * x);
}



/**
 * @brief Multiply matrix A by vector b
 * 
 * @param[in] N - size of matrix A
 * @param[in] A - matrix A with coefficients for finite difference method
 * @param[in] x - vector b with the right side of the system of equations
 * @param[out] res - vector res with result of multiplication
 */
void multiply
(
    const int N,
    const vector<double>& A, 
    const vector<double>& x, 
    vector<double>& res
) 
{
    for(int i = 0; i < N ; i++)
    {
        res[i] = 0.0;
        for(int j = 0; j < N ; j++)
        {
            res[i] += A[i * N + j] * x[j];
        }
    } 
}


double ErNorm
(
    const vector<double>& A, 
    const vector<double>& b,  
    vector<double>& x, 
    int N, 
    vector<double>& prod
)
{
    double ans = 0.0;
    double h = 1.0 / ( double(N) - 0.5 );
    double h_sq = h*h;

    multiply(N, A, x, prod);
    for (int i = 0; i < N; i++)
    {
        ans += (b[i] - prod[i]) * (b[i] - prod[i]) * h ;
    }

    

    return sqrt(ans) ;
}





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
    double h = 1.0 / ( double(N) - 0.5 );
    double h_sq = h*h;

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if     ( i == 0     && j == 0    ) A[i*N + j] = p +( 1.0 / h_sq);
            else if( i == 0     && j == 1    ) A[i*N + j] =    (-1.0 / h_sq);
            else if( i == 1     && j == 0    ) A[i*N + j] =    (-1.0 / h_sq);
            else if( i == N - 1 && j == N - 1) A[i*N + j] = p +( 2.0 / h_sq);
            else if( i == N - 1 && j == N - 2) A[i*N + j] =   -( 2.0 / h_sq);
            else if( i == N - 2 && j == N - 1) A[i*N + j] =    (-1.0 / h_sq);
            else if(i == j && i > 0 && i < N - 1 && j > 0 && j < N - 1) A[i*N + j] = p + (2.0 / h_sq);
            else if((std::abs(i - j) == 1) && i > 0 && i < N - 1 && j > 0 && j < N - 1) A[i*N + j] = (-1.0 /h_sq);
            else A[i*N + j] = 0.0;
        }
    }
}



/**
 * @brief Fill vector b with the right side of the system of equations
 * 
 * @param[in] N - size of matrix A
 * @param[out] b - vector b with the right side of the system of equations
 */
void fillf
( 
    const int N, 
    vector<double>& f
)
{
    
    double h = 1.0 / ( (double)(N) - 0.5 );
    double x_k = -h / 2.0;
    
    // f[0] = u(x_k);
    for (int k = 0; k < N  ; k++)

    {
        x_k += h;
        f[k] = u(x_k);
    }
    
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
    double h = 1.0 / (double(N) - 0.5 );
    double h_sq = h*h;
    // double l = p - 2.0 / (double(h)*double(h)) * (cos(M_PI * (2*m - 1) / (double)(2*N - 1)) - 1);
    // double l = p - 2.0 / (double(h)*double(h)) * (cos(PI * (2*m - 1) / (double)(2*N - 1)) - 1);
    double l = p - 2.0 /h_sq * (cos(2.0*M_PI * m / (2.0*(double)(N) - 1.0)) - 1.0);

    return l;
}











/**
 * @brief Solves the system of linear equations Ax = b using Richardson's method
 * 
 * @param[inout] x - the solution vector
 * @param[in] A - matrix A with coefficients for finite difference method
 * @param[in] b - vector b with the right side of the system of equations
 * @param[in] tau - parameter of the Richardson's method
 * @param[in] N - size of the matrix A
 * @param[in] mIter - number of iterations
 * @param[in] prod - vector prod with the result of the multiplication A * x
 */
double Richardson
(
    vector<double>& x, 
    const vector<double>& A, 
    const vector<double>& b, 
    const double tau, 
    const int N, 
    const int mIter,
    vector<double>& prod
)
{
    double norma = 0.0;

    for(int i = 0; i < N ; i++)
    {
        x[i]    = 0.0;
        prod[i] = 0.0;
    }

    for(int m = 0; m < mIter; m++)
    {
        multiply(N, A, x, prod);
        for(int i = 0; i < N ; i++)
        {
            x[i] = x[i] - tau * prod[i] + tau * b[i];
        }
    }
    
    return ErNorm(A, b, x, N, prod);

}











/**
 * @brief Print a portion of the matrix and a vector
 * 
 * This function prints the elements of a matrix and a vector up to a certain limit.
 * If the number of columns 'n' is less than or equal to 10, it prints all columns.
 * Otherwise, it prints the first 5 columns. Each row of the matrix is printed
 * alongside the corresponding element from vector_b.
 * 
 * @param matrix - The matrix to be printed.
 * @param n - The number of columns in the matrix.
 * @param m - The number of rows in the matrix.
 * @param vector_b - The vector to be printed alongside the matrix rows.
 */
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