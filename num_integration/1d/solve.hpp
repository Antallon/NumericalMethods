#pragma once
#include <cmath>
#include <vector>
#include <functional>

using namespace std;

double Integral_Simpson(const double a, const double b, function<double(double)> f);
double Integral_Gauss(const double a, const double b, function<double(double)> f);
double Integral_Simpson_Composite(const int N, const double a, const double b, function<double(double)> f);
double Integral_Gauss_Composite(const int N, const double a, const double b, function<double(double)> f);


