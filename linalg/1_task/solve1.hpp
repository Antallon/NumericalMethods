#pragma once
#include <cmath>
#include <vector>
#include <functional>

using namespace std;
double u(double x);
void fillmat(const int N, const double p, vector<double>& A );
void fillf(const int N, vector<double>& f);
void multiply(const int N,const vector<double>& A, const vector<long double>& x, vector<double>& res);
void make_b_from_A(vector<double>& A, int N, vector<double>& b);

long double compute_eigen_vector(const int N, const int k, const int m );
long double compute_eigen_value (const int N, const int m, const double p);
long double compute_coeff_d (const int N, const int m, const vector<double>& f);
void compute_solution_y(const int N, const double p, const vector<double>& f, vector<long double>& y );
double Richardson(vector<double>&x, const vector<double>& A, const vector<double>& b, const double tau, const int N, const int mIter, vector<double>& prod);
void printMatrix(vector<double>& matrix,int n, int m,  vector<double>& vector_b);
// double kek( vector<double> f, int N, double p);

void fillf2(const vector<double>& A,const int N,vector<double>& f,const vector<double>& x_random);


