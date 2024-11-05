#pragma once
#include <functional>
#include <vector>
using namespace std;


double u(double x);
void computePoints(int N, vector<double> &vec, vector<double>& p);
void computeTrig(int N, double x, vector<double> &sx);
double integrate(int N, vector<double> &x, function<double(double)> f);
void coeff(int N, vector<double> &x, vector<double> &coef);
double aprox(int N, double x_k, vector<double>& x, vector<double>& coef);
double compute_zero(int N, vector<double>& x);
