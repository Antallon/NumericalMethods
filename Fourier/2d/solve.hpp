#pragma once
#include <functional>
#include <vector>
#include <iterator>
using namespace std;


double u(double x, double y);
void computePoints(int N, vector<double> &vec, vector<double>& p);
void computeTrig(int N, double x, vector<double> &cx);
void coeff(int N , vector<double>& x, vector<double>::iterator coef, vector<double>::iterator u);
// double aprox(int N, double x_k, vector<double>& x, vector<double>& coef);
double compute_zero(int N, vector<double>& u);
void Ufull(int N, vector<double>& U, vector<double>& points);
void Cfull(int N, vector<double>& x,  vector<double>& C, vector<double>& U);
void Dfull(int N, vector<double>& x,  vector<double>& C, vector<double>& D);
double aprox2 (int N, double x, double y, vector<double>& points, vector<double>& U, vector<double>& D);


