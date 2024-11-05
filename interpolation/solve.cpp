#include "solve.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;


double u
(
    const double x
)
{
    // return pow(x,10) + 5.0*pow(x,5) - 2.0*pow(x,2) + 7.0*x -15.0;
    // return x*x - 3*x + 2;
    // return cos(x);
//     return std::fabs(x);
    return 1.0/(1.0 + 25*x*x);
    
    
}


/**
 * @brief Computes the norm of a matrix.
 * 
 * The norm of a matrix is the maximum of the sums of the absolute values of the
 * elements in each row.
 * 
 * @param matrix a vector of doubles representing a matrix
 * @param n the number of rows in the matrix
 * 
 * @return the norm of the matrix
 */
double matnorm
(
    const vector<double>& matrix, 
    const int n
)
{
    // Норма матрицы - максимальная сумма среди сумм элементов по строкам
    double mem = 0.0;
    double sum;
    for (int i = 0; i < n; i++)
    {   
        sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum+= std::fabs(matrix[i * n + j]);
        }
        if(sum>mem){
            mem = sum;
        }
    }
    return mem;
}



/**
 * @brief Makes a matrix triangular.
 * 
 * This function makes a matrix triangular in place. The matrix is given as a
 * vector of doubles, and the number of rows is given as an integer. The
 * function returns 1 on success, and -1 on failure.
 * 
 * @param n the number of rows in the matrix
 * @param A a vector of doubles representing a matrix
 * @param b a vector of doubles representing the right-hand side of a linear system
 * @param norma the norm of the matrix
 */
int triangle
(
    const int n, 
    vector<double>& A, 
    vector<double>& b,
    const double norma
)
{
    // printMatrix(A,n,b);

    for (int row = 0; row < n - 1; row++)                       //Идём по строкам
    {
        

        int nrow= row*n;
        // int maxCol = findMax(n, A, row);                        //Находим максимальный элемент в этой строке
        // if (maxCol != row)
        // {
        //     swapCol(n, A, maxCol, row);                         //Меняем столбцы местами
        //     printMatrix(A,n,b);
        // }
        for (int down_row = row + 1; down_row < n; down_row++)  //Идём по строкам,которые ниже текущей
        {
            int ndown_row =down_row * n;
            if (std::fabs(A[nrow + row]) <= norma * EPS)        //Если значения в матрице очень маленькие 
            {
                return -1;
            }                                                      
            double ratio = -A[ndown_row + row] / A[nrow + row]; // Находим коэфф,с которым будем вычитать из downrow row
            for (int i = row ; i < n; i++)
            {
                A[ndown_row + i] += ratio * A[nrow + i];        // Вычитаем из нижних строк текущую
            }
            b[down_row] += ratio * b[row];                      //также меняем правую часть
        }
        // printMatrix(A,n,b);


    }
    return 1;
}




/**
 * @brief Back-substitution.
 * 
 * This function performs back-substitution of a triangular matrix. The matrix is
 * given as a vector of doubles, and the number of rows is given as an integer.
 * The function returns 1 on success, and -1 on failure.
 * 
 * @param n the number of rows in the matrix
 * @param A a vector of doubles representing a triangular matrix
 * @param b a vector of doubles representing the right-hand side of a linear system
 * @param x a vector of doubles that will contain the solution of the linear system
 * @param norma the norm of the matrix
 */
int reverse
(
    const int n, 
    const vector<double>& A, 
    const vector<double>& b, 
          vector<double>& x,
    const double norma
)
{ // Обратный ход метода Гаусса
    
    for (int i = n - 1; i >= 0; i--)
    {

        if (fabs(A[i*n + i]) <= norma * EPS)
        {
            return -1;
        }

        x[i] = b[i];

        for (int j = i + 1; j < n; j++)
        {
            x[i] -= A[i*n + j] * x[j]; // из b_i вычитаем те что нашли и делим на элемент
        }
        x[i] /= A[i*n + i];
        
    }
    
    return 1;
}


/**
 * @brief Prints a matrix and a vector in a convenient format.
 * 
 * This function prints a matrix and a vector in a format that is easy to read.
 * The matrix is given as a vector of doubles, and the number of rows is given as
 * an integer. The vector is given as a vector of doubles. The function prints the
 * matrix and vector to the standard output.
 * 
 * @param matrix a vector of doubles representing a matrix
 * @param n the number of rows in the matrix
 * @param vector_b a vector of doubles representing a vector
 */
void printMatrix
(
    vector<double>& matrix, 
    int n , 
    vector<double>& vector_b
)
{
        int l;
        int m = n;
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
 * @brief Solves a system of linear equations using Gaussian elimination.
 * 
 * This function solves a system of linear equations using Gaussian elimination.
 * The function takes as input a matrix A, a vector y, and a vector coeff. The
 * matrix A is given as a vector of doubles, and the number of rows is given as
 * an integer. The vector y is given as a vector of doubles. The function writes
 * the solution into the vector coeff.
 * 
 * @param n the number of rows in the matrix A
 * @param A a vector of doubles representing a matrix
 * @param y a vector of doubles representing a vector
 * @param coeff a vector of doubles representing the solution
 * @return 0 if the function was successful, -1 otherwise
 */
int solve
(
    const int n, 
    vector<double>& A, 
    vector<double>& y, 
    vector<double>& coeff
//    const vector<double>& nodes)
       )
{
    int a, c;
    double norma = matnorm(A, n);
    // cout<<"nodes:" << endl;
    // for( size_t i  = 0; i < nodes.size(); i++) cout<<nodes[i]<<endl;
    // printMatrix(A,n,y);


    a = triangle(n, A, y, norma);
    c = reverse(n, A, y, coeff, norma);
    
    // cout<<"coeff:" << endl;
    // for( size_t i  = 0; i < coeff.size(); i++) cout<<coeff[i]<<endl;

   
    if (a == -1 || c == -1)
    {
        return -1;
    }
    // printMatrix(A,n,y);
    return 0;
}



