#include "solve.hpp"
#include "polynom.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;

double L
(
    const int n, 
    const double x, 
    const vector<double>& nodes, 
    const vector<double>& y
)
{

    double sum = 0.0;
    for(int i=0; i < n; i++)
    {
        double mult = 1.0;
        for(int j=0; j<n; j++)
        {
            if(j==i) continue;
            else
            {
                mult *= (x - nodes[j])/(nodes[i] - nodes[j]);
            }
            

        }
        sum += mult * y[i];

    }
    return sum;
}

double P
(
    const int n, 
    const double x, 
    const vector<double>& coeff
)
{
    double sum = 0.0;
    for(int i = 0; i < n; i++)
    {  
        sum += coeff[i] * pow(x,i);
    }
    return sum;    

}



void Pmatfull
(
    const int n, 
    vector<double>& A, 
    const vector<double>& nodes
)
{
    for(int i=0; i<n; i++)
    {
        double p = 1.0;
        for(int j=0; j<n; j++)
        {
            
            A[i*n + j] = p;
            p *= nodes[i];
        }
    }
}