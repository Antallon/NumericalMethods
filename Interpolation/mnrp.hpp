#pragma once
#include <cmath>
#include <vector>
using namespace std;




double mnrp(const int n, const double x, vector<double>& coeff);
void ImproveVallePuusenApproximation(const int n,const int a,const int b, vector<double>& B, vector<double>& nodes, vector<double>& coeff,  vector<double>& points);
int VallePuusen(const int n,const int a,const int b, vector<double>& B, vector<double>& nodes, vector<double>& coeff,  vector<double>& points);
void MNRPfull(const int n, vector<double>& B, const vector<double>& nodes);
