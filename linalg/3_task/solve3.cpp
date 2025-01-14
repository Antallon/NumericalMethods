#include "solve3.hpp"
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
    return cos(2 * M_PI * x);
    // return (exp(x) - exp(1) * x - 1) * (log(1 + pow(x, 2)) - 2 * log(1 + x)) * (x * x - 2 * x);

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
    const vector<double>& x, 
    vector<double>& res
) 
{
    res[0] = 0.0;
    for(int i = 0; i < N ; i++)
    {
        res[i+1] = 0.0;
        for(int j = 0; j < N ; j++)
        {
            res[i+1] += A[i * (N) + j] * x[j+1];
        }
    } 
}





/**
 * @brief Compute the norm of the difference between the exact and numerical solutions
 * 
 * @param[in] A - matrix A with coefficients for finite difference method
 * @param[in] b - vector b with the right side of the system of equations
 * @param[in] x - vector x with the numerical solution
 * @param[in] N - size of matrix A
 * @param[in] prod - vector prod with the result of multiplication A*x
 * 
 * @return norm of the difference between the exact and numerical solutions
 */
double ErNormInf(const vector<double>& A, const vector<double>& b, vector<double>& x, int N, vector<double>& prod)
{
    double ans = 0;

    multiply(N, A, x, prod); //в prod- вектор результата умножения, первый элемент мусор дальше N чисел результат умножения

    for (int i = 0; i < N; i++)
    {
        if(fabs((b[i+1] - prod[i+1])) > ans) // в b все так же нулевой элемент мусор
        { 
            ans = fabs((b[i+1] - prod[i+1]));
        }
    }

    return ans;
}




/**
 * @brief Fill matrix A with coefficients for finite difference method
 * 
 * @param[in] N - size of matrix A
 * @param[in] p - parameter for matrix A filling
 * @param[in,out] A - matrix A with coefficients for finite difference method
 * @param[in] flag - flag for filling matrix A
 */

void fillMat
(
    int N,
    double p, 
    vector<double>& A, 
    bool flag
) 
{
    double h = 1.0 / (double)(N - 0.5);
    double h_sq = h*h;
    p=0.0;
    for(int i = 0; i < N; i++)
    {
        p+= 1.0 + (double)sin((double)M_PI*(i+1)*h)*(double)sin((double)M_PI*(i+1)*h);
    }
    p *=h;



    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if(flag == true){
                p = 1.0 + (double)sin((double)M_PI*(j+1)*h)*(double)sin((double)M_PI*(j+1)*h);
                // p = 1.0;

                // h_sq = M_PI * M_PI / ((double)(N - 0.5) * (double)(N - 0.5));
            }
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
 * @brief Fill the vector f with values of function u(x) on the uniform grid
 * 
 * @param[in] N - size of the problem
 * @param[out] f - vector of function values
 */
void fillf
( 
    const int N, 
    vector<double>& f
)
{
    double h = 1.0 / ( (double)(N) - 0.5 );
    double x_k = -h / 2.0;
    f[0] = u(x_k);
    for (int k = 1; k < N  ; k++)
    {
        x_k += h;
        f[k] = u(x_k);
    }
    f[N] = u(1.0);
}

void fillf2
(
    const vector<double>& A,
    const int N,
    vector<double>& f,
    const vector<double>& x_random
)
{
    for(int i = 0; i < N ; i++)
    {
        f[i+1] = 0.0;
        for(int j = 0; j < N ; j++)
        {
            f[i+1] += A[i * (N) + j] * x_random[j+1];
        }
    } 
}


























/**
 * @brief Compute eigen vector of the matrix A
 * 
 * @param[in] N - size of matrix A
 * @param[in] k - number of eigen vector
 * @param[in] m - number of element in eigen vector
 * @return eigen vector of the matrix A
 */
double compute_eigen_vector
(
    const int N, 
    const int k, 
    const int m 
)
{
    return cos(M_PI * m*(2.0*k - 1.0) / (2.0*(double)N - 1.0));
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
    const vector<double>& f
)
{
    double h = 1.0 / ((double)N - 0.5);
    double sum = 0.0;
    double sum0 = 0.0;
    // for (int k = 1; k < N; k++)
    // {
    //     sum0 +=  compute_eigen_vector(N, k, m) * compute_eigen_vector(N, k, m) * (h);
    // }
    // sum0 +=  compute_eigen_vector(N, N, m) * compute_eigen_vector(N, N, m) * (h/2.0);
    // cout<<"sum0 = " << sum0<<endl;
    sum0 = 2.0;
    for (int k = 1; k < N; k++)
    {
        sum +=   f[k] * compute_eigen_vector(N, k, m) * (h);
    }
    sum +=   f[N] * compute_eigen_vector(N,N,m) * (h / 2.0);
    // cout<<"sum"<<sum<<endl;
    return sum;
}



/**
 * @brief Compute the solution vector y for the problem
 * 
 * @param[in] N - size of the problem
 * @param[in] p - parameter of the problem
 * @param[out] f - vector of function values
 * @param[in] y -  values solution vector
 */
void compute_solution_y
(
    const int N, 
    const double p,
    const vector<double>& f,
    vector<double>& y,
    vector<double>& prod2,
    const vector<double>& B

)
{
    for(int k = 0; k < N + 1 ; k++)
    {
        y[k] = 0.0;
        for(int m = 1; m < N; m++)
        {
            y[k] += 2.0*compute_coeff_d(N, m, f) / compute_eigen_value(N, m, p) * compute_eigen_vector(N, k, m) ;
        }
        y[k] += compute_coeff_d(N, N, f) / compute_eigen_value(N, N, p) * compute_eigen_vector(N, k, N) ;
        y[k] += compute_coeff_d(N, 0, f) / compute_eigen_value(N, 0, p) * compute_eigen_vector(N, k, 0) ;

    }
    multiply(N, B, y, prod2);
    for(int i = 1; i <=N; i++)
    {
       
    }
    
}




































//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







/**
 * @brief Find the q value for the method of simple iteration
 * 
 * @param[in] A - matrix of coefficients
 * @param[in] N - size of the matrix
 * 
 * @return q value
 */
double Findq
(
    const vector<double>& A, 
    int N
)
{
    //из диагонального преобладания
    double sum = 0;
    double max_q = 0;
    double q_temp;
    int i, j;


    for (i = 0; i < N; i++){
        sum = 0.0;
        for (j = 0; j < N; j++){
            if(j != i){
                sum += fabs(A[i * N + j]);
            }
        }

        q_temp = sum / fabs(A[i * N + i]);
        if(q_temp > 1){
            printf("bad q, would haven't conv\n");
        }
        if(i == 0){
            max_q = q_temp;
        }
        if (q_temp > max_q){
            max_q = q_temp;
        }
    }
    return max_q;
}








// /**
//  * @brief Method of simple iteration
//  * 
//  * @param[inout] x - values of unknown vector
//  * @param[in] A - matrix of coefficients
//  * @param[in] B - matrix of coefficients
//  * @param[in] b - vector of right-hand side values
//  * @param[in] tau - step size of simple iteration method
//  * @param[in] N - size of the problem
//  * @param[in] mIter - number of iterations
//  * @param[inout] prod - temporary vector for multiplication
//  * @param[inout] prod1 - temporary vector for storing y_k
//  * @param[in] p - parameter of the problem
//  * @return the norm of the residual at the last iteration
// */
// double BSolver
// (
//     vector<double>& x, 
//     const vector<double>& A, 
//     const vector<double>& B, 
//     const vector<double>& b, 
//     double tau, 
//     int N, 
//     int mIter, 
//     vector<double>& prod, 
//     vector<double>& prod1,
//     vector<double>& prod2,
//     double p
// )
// {
//     for (int k = 0; k < N+1; k++)
//     {
//         x[k]     = 0.0;
//         prod[k]  = 0.0;
//         prod1[k] = 0.0;
//         prod2[k] = 0.0;
//     }
//     for (int m = 0; m < mIter; m++)
//     {
//         multiply(N, A, x, prod);
//         // cout<<"b- AX for iter m = "<<m<<endl;
//         for (int j = 1; j < N+1; j++)
//         {
//             prod[j] = b[j] - prod[j];
//             // cout<<prod[j]<<endl;
//         }
//         // cout<<endl;
//         compute_solution_y(N, p, prod, prod1,prod2, B); //prod1 - y_k из системы By_k = b - Ax_k
//         for (int i = 1; i <= N; i++)
//         {
//             x[i] = x[i] + tau * prod1[i];
//         }
//     }
//     return ErNormInf(A, b, x, N, prod);
// }





/**
 * @brief Method of simple iteration
 * 
 * @param[inout] x - values of unknown vector
 * @param[in] A - matrix of coefficients
 * @param[in] B - matrix of coefficients
 * @param[in] b - vector of right-hand side values
 * @param[in] tau - step size of simple iteration method
 * @param[in] N - size of the problem
 * @param[in] mIter - number of iterations
 * @param[inout] prod - temporary vector for multiplication
 * @param[inout] prod1 - temporary vector for storing y_k
 * @param[in] p - parameter of the problem
 * @return the norm of the residual at the last iteration
*/
double BSolver
(
    vector<double>& x, 
    const vector<double>& A, 
    const vector<double>& B, 
    const vector<double>& b, 
    double tau, 
    int N, 
    int mIter, 
    vector<double>& prod, 
    vector<double>& prod1,
    vector<double>& prod2,
    double p
)
{
    for (int k = 0; k < N+1; k++)
    {
        // x[k]     = 0.0;
        prod[k]  = 0.0;
        prod1[k] = 0.0;
        prod2[k] = 0.0;
    }
    
    multiply(N, A, x, prod);
    for (int j = 1; j < N+1; j++)
    {
        prod[j] = b[j] - prod[j];
    }
        

    compute_solution_y(N, p, prod, prod1,prod2, B); //prod1 - y_k из системы By_k = b - Ax_k
    for (int i = 1; i <= N; i++)
    {
        x[i] = x[i] + tau * prod1[i];
    }
    return ErNormInf(A, b, x, N, prod);
}









