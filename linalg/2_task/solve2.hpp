#pragma once
#include <cmath>
#include <vector>
#include <functional>

using namespace std;
double u(double x);
void multiply(const int N,const vector<double>& A, const vector<double>& x, vector<double>& res);
double ErNorm(const vector<double>& A, const vector<double>& b,  vector<double>& x, int N, vector<double>& prod);
void fillmat(const int N, const double p, vector<double>& A );
void fillf(const int N, vector<double>& f);
double Richardson(vector<double>&x, const vector<double>& A, const vector<double>& b, const double tau, const int N, const int mIter, vector<double>& prod);
double compute_eigen_value (const int N, const int m, const double p);
void printMatrix(vector<double>& matrix,int n, int m,  vector<double>& vector_b);

