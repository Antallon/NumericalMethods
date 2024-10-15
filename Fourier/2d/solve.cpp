#include <math.h>
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
using namespace std;



double u(double x, double y)
{
    // return cos(2 * M_PI * x)*cos(2 * M_PI * y);
    // return sin(2 * M_PI * x)*sin(2 * M_PI * y);
    return cos(2 * M_PI * x)*sin(2 * M_PI * y);
    // return cos(2 * M_PI * x)*cos(2 * M_PI * y) * cos(2 * M_PI * x)*cos(2 * M_PI * y);


    // return x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y);
    // return x*x*(1-x*x) * sin(M_PI * x)*y*y*(1-y*y) * sin(M_PI * y);
    // return x*y;
}


void computePoints(int N, vector<double>& vec,vector<double>& p)
{
    double h = 1.0 / ( N - 0.5 );
    double x_k = -h / 2;
    vec[0] = x_k;
    p[0] = x_k;

    for(int i=1; i < N; i++)
    {
        x_k += h;
        vec[i] = x_k;
    }

    for(int j = 1; j < (N-1)*3 + 1; j++)
    {
        p[j] = p[j-1] + h/3;
    }
}


void computeTrig(int N, double x, vector<double>& cx)
{
    for(int m=0; m < N; m++)
    {
        cx[m] = cos(M_PI * m * x);
    }
}



double compute_zero(int N, vector<double>& u)
{
    double h = 1.0 / ( N - 0.5 );
    double zero = 0.0;
    for (int i = 0; i < N; i++)
    {
        for(int j=0; j < N; j++)
        {
            zero += u[i * N +j] * h * h ;
        }
    }

    return zero; 
}


void coeff(int N , vector<double>& x, vector<double>::iterator coef, vector<double>::iterator u)
{
    double h = 1.0 / ( N - 0.5 );
    
    for(int i=0; i < N; i++)
    {
        *(coef + i) = 0.0;
        for(int j=1; j < N; j++)
        {
            *(coef + i) += *(u+j) * h * cos(M_PI * i * x[j]);
        }
        if(i!=0) *(coef + i) *= 2;
    }

}

void Ufull(int N, vector<double>& U, vector<double>& points)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            
            U[i * N + j] =  u(points[i], points[j]);
            
        }
    }
}


void Cfull(int N, vector<double>& x,  vector<double>& C, vector<double>& U)
{
  

    for(int i = 0; i < N; i++)
    {
        // C[i*N] = 0;
        coeff(N,x,C.begin() + i * N, U.begin() + i * N);
    }
}

void Dfull(int N, vector<double>& x,  vector<double>& C, vector<double>& D)
{
    for(int i = 0; i < N; i++)
    {
        vector <double> col(N);
        // D[i*N] = 0;
        for(int h = 0; h < N; h++)
        {
            col[h] = C[h * N + i];
        }
        coeff(N, x, D.begin() + i * N, col.begin());
    }
}


double aprox2 (int N, double x, double y, vector<double>& points, vector<double>& U, vector<double>& D)
{
    vector <double> cx(N);
    vector <double> cy(N);
    computeTrig(N, x, cx);
    computeTrig(N, y, cy);
    // double z = compute_zero(N, U);
    // cout<<z<<endl;
    double z = 0.0;
    for(int i=0; i < N; i++)
    {
        for(int j=0; j < N; j++)
        {
            z += D[i * N + j] * cx[i] * cy[j];
        }
    }
    return z;
}








