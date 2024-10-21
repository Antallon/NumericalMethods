#include <cmath>
#include <vector>
using namespace std;

double u(const double x);
double P(const int n, const double x, const vector<double>& coeff);
double L(const int n, const double x, const vector<double>& nodes, const vector<double>& y);
double mnrp(const int n, const double x, vector<double>& coeff);


void Pmatfull(const int n, vector<double>& A, const vector<double>& nodes);
void MNRPfull(const int n, vector<double>& B, const vector<double>& nodes);

double matnorm(const vector<double>& matrix, const int n);

// void swapCol(const int n, vector<double>& A, const int col1, const int col2);
// int findMax(const int n, const vector<double>& A, const int row);
int triangle(const int n, vector<double>& A, vector<double>& b,const double norma);
int reverse(const int n, const vector<double>& A, const vector<double>& b, vector<double>& x,const double norma);
int solve(const int n, vector<double>& A, vector<double>& y, vector<double>& coeff, const vector<double>& nodes);

void printMatrix(vector<double>& matrix, int n , vector<double>& vector_b);


