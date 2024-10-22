#pragma once
#include <cmath>
#include <vector>
using namespace std;



double P(const int n, const double x, const vector<double>& coeff);
double L(const int n, const double x, const vector<double>& nodes, const vector<double>& y);
void Pmatfull(const int n, vector<double>& A, const vector<double>& nodes);
