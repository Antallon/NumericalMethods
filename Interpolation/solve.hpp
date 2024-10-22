#pragma once
#include <cmath>
#include <vector>
using namespace std;

double u(const double x);
double matnorm(const vector<double>& matrix, const int n);
int triangle(const int n, vector<double>& A, vector<double>& b,const double norma);
int reverse(const int n, const vector<double>& A, const vector<double>& b, vector<double>& x,const double norma);
int solve(const int n, vector<double>& A, vector<double>& y, vector<double>& coeff, const vector<double>& nodes);
void printMatrix(vector<double>& matrix, int n , vector<double>& vector_b);



// void swapCol(const int n, vector<double>& A, const int col1, const int col2)
// {
//     // Функция меняющая колонки местами
//     for (int i = 0; i < n; i++)
//     {   
//         std::swap(A[i*n + col1], A[i*n + col2]);
//     }
// }

// int findMax(const int n, const vector<double>& A, const int row)
// {
//     // Функция получает на вход строку, находит в ней максимальный элемент
//     //и возвращает столбец,в котором стоит максимальный элемент
//     int maxCol = row;
//     for (int i = row + 1; i < n; i++)
//     { 
//     // идем вправо по колонкам и ищем максимальный элемент в строке
//         int row_n = row * n;
//         if (fabs(A[row_n + i]) > fabs(A[row_n + maxCol]))
//         {
//             maxCol = i;
//         }
//     }
//     return maxCol;
// }