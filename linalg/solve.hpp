#pragma once
#include <cmath>
#include <vector>
#include <functional>

using namespace std;
void fillmat(const int N, const double p, vector<double>& A );
void fillf(const vector<double>& A, const int N, vector<double>& f);
void multiply(const int N,const vector<double>& A, const vector<double>& b, vector<double>& res);
double vecnorm(const vector<double>& v);
void vecdiff(const vector<double>& v1, const vector<double>& v2, vector<double>& diff);

double compute_eigen_vector(const int N, const int m, const int k );
double compute_eigen_value (const int N, const int m, const double p);
double compute_coeff_d (const int N, const int m, const double p, const vector<double>& f);
void compute_solution_y(const int N, const double p, const vector<double>& f, vector<double>& y );
void Richardson(vector<double>&x, const vector<double>& A, const vector<double>& b, const double tau, const int N, const int mIter, vector<double>& prod);
void printMatrix(vector<double>& matrix,int n, int m,  vector<double>& vector_b);

