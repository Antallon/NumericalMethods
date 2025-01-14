#pragma once
#include <cmath>
#include <vector>
#include <functional>

using namespace std;
double u(double x);

double ErNormInf(const vector<double>& A, const vector<double>& b, vector<double>& x, int N, vector<double>& prod);
void multiply(const int N,const vector<double>& A, const vector<double>& b, vector<double>& res);
void fillMat(int N,double p, vector<double>& A, bool flag);
void fillf(const int N, vector<double>& f);

double compute_eigen_vector(const int N, const int m, const int k );
double compute_eigen_value (const int N, const int m, const double p);
double compute_coeff_d (const int N, const int m, const double p, const vector<double>& f);
void compute_solution_y(const int N, const double p, const vector<double>& f, vector<double>& y, vector<double>& prod2, const vector<double>& B);
void printMatrix(vector<double>& matrix,int n, int m,  vector<double>& vector_b);

double BSolver(vector<double>& x, const vector<double>& A, const vector<double>& B, const vector<double>& b, double tau, int N, int mIter, vector<double>& prod, vector<double>& prod1,vector<double>& prod2, double p);
double Findq(const vector<double>& A, int n);


void fillf2(const vector<double>& A,const int N,vector<double>& f,const vector<double>& x_random);
