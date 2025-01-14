#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;


double u(double x);
double div_diff(const vector<int>& basis, const vector<double>& net, const vector<double>& y, int start, int end);
double div_diff_phi(const vector<int>& basis, const vector<double>& net, const vector<double>& phi, int start, int end);
double compute_h(const vector<int>& basis, const vector<double>& net, const vector<double>& fun, const vector<double>& phi, int n);
double compute_polynom_in_point(const vector<double>& fun, const vector<int>& basis, const vector<double>& net, const vector<double>& phi, double h, double x, int n);
void compute_polynom_values(vector<double>& poly, const vector<double>& fun, const vector<int>&basis, const vector<double>& net, const vector<double>& phi, int n, int N);
double compute_h_sigma(const vector<int>& basis, const vector<double>& fun, const vector<double>& poly, int n);
double compute_phi_sigma(const vector<double>& fun, const vector<double>& poly, int N);
void transform_basis(vector<int>& basis, const vector<double>& fun, const vector<double>&poly, double phi_sigma, int n, int N);
void vallee_pussen_method(int n, int N,vector<int>& basis,vector<double>& phi, vector<double>& poly,const vector<double>& nodes,const vector<double>& fun);


