#include <math.h>
#include <cmath>
#include <functional>
#include <vector>
using namespace std;

double u(double x)
{
    // return cos(2* M_PI * x)/2 + cos(M_PI * x);
    // return x*x*(1-x*x) * sin(M_PI * x); 
    // return cos(M_PI * x) * cos(M_PI * x); 
    return cos(2*M_PI * x) * sin(2*M_PI * x); 

    // return x*x*(1-x)*(1-x) ;
    // + cos(M_PI * x);
    // return cos(2 * M_PI * x);
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


void computeTrig(int N, double x, vector<double>& sx)
{
    // sx[0] = 0;
    for(int m=0; m < N; m++)
    {
        // sx[m] = sin(M_PI * m * x);
        sx[m] = cos(M_PI * m * x);
        // sx[m] = cos(2 * M_PI * m * x / (2*N-1));
    }
}


double integrate(int N, vector<double>& x, function<double(double)> f)
{
    double h = 1.0 / ( N - 0.5 );
    double integral = 0.0;
    for(int i=1; i < N; i++)
    {
        integral+= f(x[i]) * h;
    }
    return integral;
}

double compute_zero(int N, vector<double>& x)
{
    double zero =  integrate(N ,x , [](double a){
        return u(a);
    });
    return zero;
}


void coeff(int N , vector<double>& x, vector<double>& coef)
{

    for(int i=1; i < N; i++)
    {
        coef[i] = 2 * integrate(N, x, [i](double a){
            return u(a) * cos(M_PI * i * a);
        });
    }

}

double aprox(int N, double x_k, vector<double>& x, vector<double>& coef)
{
    vector<double> sx(N);
    computeTrig(N, x_k, sx);
    double y_k = compute_zero(N,x);
    // double y_k = 0.0;
    for(int i=1; i < N; i++)
    {
        y_k+=coef[i] * sx[i];
    }
    return y_k;
}








