#include "solve1.hpp"
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
    return x*x*(1-x)*(1-x) ;
    // + cos(M_PI * x);
    // return cos(2 * M_PI * x);
    // return cos(10*M_PI*x);
    // return sin(M_PI*x);
    // return x*x*x*x + 3*x*x*x - 1;

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
    const vector<long double>& x, 
    vector<double>& res
) 
{
    for(int i = 0; i < N ; i++)
    {
        res[i] = 0.0;
        for(int j = 0; j < N ; j++)
        {
            res[i] += A[i * (N) + j] * x[j+1];
        }
    } 
}

// void make_b_from_A(vector<double>& A, int N, vector<double>& b)
// {
//     double h = 1.0 / ( double(N) - 0.5 );
//     for (int k = 1; k < N; k++)
//     {
//         b[k] = 0.;
//         for (int j = 0; j < N+1; j++)
//         {
//             b[k] += A[k * (N+1) + j];
//         }
//     }
//     // b[0] = b[1];
//     b[0] = 0.0;
//     b[N] = 0.0;
// }





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
    // double rofl = (2 - p*p*h*h*h*h) / (2 - p*h*h);
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if     ( i == 0     && j == 0    ) A[i*N + j] = p +( 1.0 / h_sq);
            else if( i == 0     && j == 1    ) A[i*N + j] =    (-1.0 / h_sq);
            else if( i == 1     && j == 0    ) A[i*N + j] =    (-1.0 / h_sq);
         // else if( i == N - 1 && j == N - 1) A[i*N + j] = p -( 2.0 / h_sq);
            else if( i == N - 1 && j == N - 1) A[i*N + j] = p +( 2.0 / h_sq);
            // else if( i == N - 1 && j == N - 2) A[i*N + j] =    ( 2.0 / h_sq);
            else if( i == N - 1 && j == N - 2) A[i*N + j] =    -( 1.0 / h_sq); //странное место
 
            else if( i == N - 2 && j == N - 1) A[i*N + j] =    (-1.0 / h_sq);
            else if(i == j && i > 0 && i < N - 1 && j > 0 && j < N - 1) A[i*N + j] = p + (2.0 / h_sq);
            else if((std::abs(i - j) == 1) && i > 0 && i < N - 1 && j > 0 && j < N - 1) A[i*N + j] = (-1.0 /h_sq);
            else A[i*N + j] = 0.0;

            // if(i == 0 && j == 0) A[i*N + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if((i == 0 && j == 1) || (i==N-2 && j==N-1)) A[i*N + j] = -1.0 / ((double)(h) * (double)(h));
            // else if(i == N-1 && j == N - 2) A[i*N + j] = -2.0 / ((double)(h) * (double)(h));
            // else if(i == N-1 && j == N - 1) A[i*N + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if(i == 1 && j == 0)A[i*N + j] = -1.0 / ((double)(h) * (double)(h));
            // else if(i == j && i > 0 && i < N - 1 && j > 0 && j < N - 1) A[i*N + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if((std::abs(i - j) == 1) && i > 0 && i < N - 1 && j > 0 && j < N - 1) A[i*N + j] = -1.0 / ((double)(h) * (double)(h));
            // else A[i*N + j] = 0.0;
            /////////////////////////////////////////
            // if(i == 0 && j == 0) A[i*N + j] = -1.0 / h;
            // else if((i == 0 && j == 1)) A[i*N + j] = 1.0 / h;
            // if( (i==N-2 && j==N-1)) A[i*(N+1) + j] = -1.0 / ((double)(h) * (double)(h));
            // else if(i == 1 && j == 1)A[i*(N+1) + j] = p + 1.0 / ((double)(h) * (double)(h));

            // else if(i == N-1 && j == N - 2) A[i*(N+1) + j] = -2.0 / ((double)(h) * (double)(h));
            // else if(i == N-1 && j == N - 1) A[i*(N+1) + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if(i == N && j == N ) A[i*(N+1) + j] = (-2.0 / ((double)(h) * (double)(h))) + p;
            // else if(i == N && j == N - 1 ) A[i*(N+1) + j] =  2.0 / ((double)(h) * (double)(h));
            // else if(i == N - 1 && j == N  ) A[i*(N+1) + j] = - 1.0 / ((double)(h) * (double)(h));
            // else if(i == 1 && j == 0)A[i*(N+1) + j] = -1.0 / ((double)(h) * (double)(h));
            // else if(i == j && i > 0 && i < N  && j > 0 && j < N ) A[i*(N+1) + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if((std::abs(i - j) == 1) && i > 0 && i < N  && j > 0 && j < N ) A[i*(N+1) + j] = -1.0 / ((double)(h) * (double)(h));
            // else A[i*(N+1) + j] = 0.0;
            //////////////////////////////////////////




///////////////////////
            // if     (i == 0   && j == 0    ) A[i*(N+1) + j] = -1.0;// /  double(h);
            // else if(i == 0   && j == 1    ) A[i*(N+1) + j] =  1.0;// /  double(h);
            // else if(i == N-1 && j == N    ) A[i*(N+1) + j] = -1.0 / ((double)(h) * (double)(h));
            // else if(i == N   && j == N - 1) A[i*(N+1) + j] =  2.0 / ((double)(h) * (double)(h));
            // else if(i == N   && j == N    ) A[i*(N+1) + j] = -2.0 / ((double)(h) * (double)(h)) + p;
            // // else if(i == N   && j == N - 1) A[i*(N+1) + j] =  -2.0 / ((double)(h) * (double)(h));
            // // else if(i == N   && j == N    ) A[i*(N+1) + j] =  2.0 / ((double)(h) * (double)(h)) + p;
            // else if(i == 1   && j == 0    ) A[i*(N+1) + j] = -1.0 / ((double)(h) * (double)(h));
            // else if(i == j   && i > 0 && i < N  && j > 0 && j < N ) A[i*(N+1) + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if((std::abs(i - j) == 1) && i > 0 && i < N  && j > 0 && j < N ) A[i*(N+1) + j] = -1.0 / ((double)(h) * (double)(h));
            // else A[i*(N+1) + j] = 0.0;
//////////////////////////



/////////////////////////////////////////////
            // if     (i == 1   && j == 1    ) A[i*(N+1) + j] =  1.0 / ((double)(h) * (double)(h)) + p;// /  double(h);
            // else if(i == N-1   && j == N-1    ) A[i*(N+1) + j] =  rofl;
            // else if(i == j   && i > 0 && i < N  && j > 0 && j < N ) A[i*(N+1) + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if((std::abs(i - j) == 1) && i > 0 && i < N  && j > 0 && j < N ) A[i*(N+1) + j] = -1.0 / ((double)(h) * (double)(h));
            // else A[i*(N+1) + j] = 0.0;
///////////////////////////////////

            // else if(i == j   && i > 0 && i < N  && j > 0 && j < N ) A[i*(N+1) + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if((std::abs(i - j) == 1) && i > 0 && i < N  && j > 0 && j < N ) A[i*(N+1) + j] = -1.0 / ((double)(h) * (double)(h));
            // else A[i*(N+1) + j] = 0.0;


            // if(i==1 && j == 1) A[i * (N + 1) + j] = p + 1.0 / ((double)(h) * (double)(h));
            // else if(i == N - 1 && j == N-1) A[i * (N + 1) + j] = p + 2.0 / ((double)(h) * (double)(h));
            // else if(i == N - 1 && j == N) A[i * (N + 1) + j] = -2.0 / ((double)(h) * (double)(h));
            // else if(i == j && i > 0 && i < N && j > 0 && j < N) A[i * (N + 1) + j] = p + 2.0 / ((double)(h) * (double)(h));

            // else if((std::abs(i - j) == 1) && i > 0 && i < N && j > 0 && j < N) 
            // A[i * (N + 1) + j] = -1.0 / ((double)(h) * (double)(h));
            
            // else A[i * (N + 1) + j] = 0.0;
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
    
    f[0] = u(x_k);
    for (int k = 1; k < N  ; k++)
    {
        x_k += h;
        f[k] = u(x_k);
    }
    // f[N] = 0.0;
    f[N] = u(1.0);
    // for (int k = 0; k < N + 1 ; k++)
    // {
    //     f[k] = u(0.0);
    //     // if(k == 5) f[k] = 1.0;
    // }
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
long double compute_eigen_vector
(
    const int N, 
    const int k, 
    const int m 
)
{
    // dmake
    // long double ev = 0.0;
    long double ev = cos(M_PI * m*(2.0*k - 1.0) / (2.0*N - 1.0));
    // long double ev = sin(PI * k * (2 * m - 1) / (double)(2 * N - 1));
    return ev;
}



/**
 * @brief Compute eigen value of the matrix A
 * 
 * @param[in] N - size of matrix A
 * @param[in] m - number of eigen value
 * @param[in] p - parameter of the problem
 * @return eigen value of the matrix A
 */
long double compute_eigen_value
(
    const int N, 
    const int m, 
    const double p
)
{
    long double l = 0.0;
    long double h = 1.0 / (double(N) - 0.5 );
    long double h_sq = h*h;
    // double l = p - 2.0 / (double(h)*double(h)) * (cos(M_PI * (2*m - 1) / (double)(2*N - 1)) - 1);
    // double l = p - 2.0 / (double(h)*double(h)) * (cos(PI * (2*m - 1) / (double)(2*N - 1)) - 1);
    l = p - 2.0 /h_sq * (cos(2.0*M_PI * m / (2.0*(double)(N) - 1.0)) - 1.0);
    return l;
}




/**
 * @brief Compute the coefficient d for the finite difference method
 * 
 * @param[in] N - size of the problem
 * @param[in] m - index of the eigen vector
 * @param[in] f - vector of function values
 * @return coefficient d based on the given eigen vector and function values
 */
long double compute_coeff_d
(
    const int N, 
    const int m, 
    const vector<double>& f
)
{
    long double h = 1.0 / ((double)N - 0.5);
    // double h = 1.0 / ((double)N);
    long double sum = 0.0;
    long double sum0 = 0.0;
    // for (int k = 1; k < N; k++)
    // {
    //     sum0 +=  compute_eigen_vector(N, k, m) * compute_eigen_vector(N, k, m) * (h);
    // }
    // sum0 +=  compute_eigen_vector(N, N, m) * compute_eigen_vector(N, N, m) * (h/2.0);
    // sum0 +=  compute_eigen_vector(N, 0, m) * compute_eigen_vector(N, 0, m) * (h/2.0);

    // cout<<"sum0 = " << sum0<<endl;
    sum0 = 2.0;
    for (int k = 1; k < N; k++)
    {    
        sum +=  double(f[k]) * double(compute_eigen_vector(N, k, m)) * double(h);        
    }
   
    sum +=   double(f[N]) * double(compute_eigen_vector(N,N,m)) * double((h / 2.0));
    // sum +=   2.0*double(f[0]) * double(compute_eigen_vector(N,0,m)) * double((h / 2.0));

   
    // cout<<"sum"<<sum<<endl;
    // cout<<"END////////"<<endl;
    // cout<<endl;
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
    vector<long double>& y

)
{
    for(int k = 0; k < N + 1 ; k++)
    {
        y[k] = 0.0;
        for(int m = 1; m < N; m++)
        {
            if(m == 0 || m == N)
            {
                // cout<<"m = "<<m<<endl;
                // cout<<"k = "<<k<<endl;
                // cout<<"compute_eigen_vector = "<<compute_eigen_vector(N, k, m)<<endl;
                // cout<<"compute_eigen_value = "<<compute_eigen_value(N, m, p)<<endl;
                // cout<<"compute_coeff_d = "<<compute_coeff_d(N, m, f)<<endl;
                // cout<<"plus = "<< compute_coeff_d(N, m, f) / compute_eigen_value(N, m, p) * compute_eigen_vector(N, k, m)<<endl ;
                // cout<<"y_k before = "<<y[k]<<endl;
            }
            y[k] += 2.0 * compute_coeff_d(N, m, f) / compute_eigen_value(N, m, p) * compute_eigen_vector(N, k, m) ;

        }
        y[k] += compute_coeff_d(N, N, f) / compute_eigen_value(N, N, p) * compute_eigen_vector(N, k, N) ;
        y[k] += compute_coeff_d(N, 0, f) / compute_eigen_value(N, 0, p) * compute_eigen_vector(N, k, 0) ;


    }


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