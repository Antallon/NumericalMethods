#include "solve.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#define EPS 1e-16
using namespace std;


double u(const double x)
{
    return pow(x,10) + 5.0*pow(x,5) - 2.0*pow(x,2) + 7.0*x -15.0;
    // return x*x - 3*x + 2;
    // return cos(x);
    // return std::fabs(x);
    // return 1.0/(1.0 + x*x);
    
    
}

double P(const int n, const double x, const vector<double>& coeff)
{
    double sum = 0.0;
    for(int i = 0; i < n; i++)
    {  
        sum += coeff[i] * pow(x,i);
    }
    return sum;    cout << endl;

}
double L(const int n, const double x, const vector<double>& nodes, const vector<double>& y)
{
    // cout<<"====================="<<endl;
    // cout<< " x = " << x << endl;
    // cout<<endl;
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
        // cout<<" Phi_"<<i<<" = "<< mult<<endl;
        // cout<<"sum before adding = "<< sum<<endl;

        sum += mult * y[i];
        // cout<<"mult * y_i = "<<mult<<" *"<< y[i]<<" = "<< mult * y[i]<<endl;
        // cout<<"sum = " << sum<<endl;
        // cout<<endl;

    }
    return sum;
}


void matfull(const int n, vector<double>& A, const vector<double>& nodes)
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
double matnorm(const vector<double>& matrix, const int n)
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



void swapCol(const int n, vector<double>& A, const int col1, const int col2)
{
    // Функция меняющая колонки местами
    for (int i = 0; i < n; i++)
    {   
        std::swap(A[i*n + col1], A[i*n + col2]);
    }
}

int findMax(const int n, const vector<double>& A, const int row)
{
    // Функция получает на вход строку, находит в ней максимальный элемент
    //и возвращает столбец,в котором стоит максимальный элемент
    int maxCol = row;
    for (int i = row + 1; i < n; i++)
    { 
    // идем вправо по колонкам и ищем максимальный элемент в строке
        int row_n = row * n;
        if (fabs(A[row_n + i]) > fabs(A[row_n + maxCol]))
        {
            maxCol = i;
        }
    }
    return maxCol;
}

int triangle(const int n, vector<double>& A, vector<double>& b,const double norma)
{
    // printMatrix(A,n,b);

    for (int row = 0; row < n - 1; row++)                       //Идём по строкам
    {
        

        int nrow= row*n;
        int maxCol = findMax(n, A, row);                        //Находим максимальный элемент в этой строке
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




int reverse(const int n, const vector<double>& A, const vector<double>& b, vector<double>& x,const double norma)
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


void printMatrix(vector<double>& matrix, int n , vector<double>& vector_b){
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



   




int solve(const int n, vector<double>& A, vector<double>& y, vector<double>& coeff, const vector<double>& nodes)
{
    int a, c;
    double norma = matnorm(A, n);
    // cout<<"nodes:" << endl;
    // for( size_t i  = 0; i < nodes.size(); i++) cout<<nodes[i]<<endl;
    matfull(n,A,nodes);
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


