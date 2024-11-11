#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>                                               //Для setw
#include "solve.hpp"
#include <cmath>

#include <math.h>
#include <cmath>
using namespace std;
#define pi 3.14159265358979323846



int main(int argc,  char **argv)
{   
    if(argc == 2) cout<<endl;
    int N, mIter;
    double p, tau,m, M, q, norma_zero;
    double mult = 1.0;
    ofstream File("result.txt");
    cout<<"Process start "<<endl<<endl;

    N     = std::stoi(argv[1]);
    mIter = std::stoi(argv[2]);
    p     = std::stod(argv[3]);
   
    if (!File)
    {
        cerr << "Error: Unable to open file" << endl;
        return -1;
    }
    





    vector<double> A((N+1)*(N+1));
    vector<double> x(N + 1);
    vector<double> y(N+1);
    vector<double> f(N+1);
    vector<double> prod(N+1);
    vector<double> diff_zero(N+1);
    vector<double> diff(N+1);

    m = compute_eigen_value(N, 1, p);
    M = compute_eigen_value(N, N, p);
    tau = 2.0 / (m + M);
    q = (M - m) / (M + m);
    cout<<"q = "<<q<<endl;
    cout<<"m = "<<m<<endl;
    cout<<"M = "<<M<<endl;


//  Task 1
    fillmat(N, p, A);
    fillf(A, N, f);
    compute_solution_y(N, p, f, y);
    multiply(N, A, y, prod);

   



    
    
    for(int k = 0; k < mIter; k++)
    {
        Richardson(x, A, f, tau, N, k, prod);

        if(k == 0)
        {
            vecdiff(y, x, diff_zero);
            norma_zero = vecnorm(diff_zero);
            for(int i = 0; i < N + 1; i++) cout<< x[i]<<" ";
            cout<<endl<<endl;
            cout<<norma_zero;
        }
        else 
        {
            vecdiff(y, x, diff);
        }
        File
        <<setw(20)<< k
        <<setw(20)<<vecnorm(diff)
        <<setw(20)<<mult * norma_zero
        << endl;
        mult *= q;
    }
    
    

















    cout<<"The solution for method Fourier is: "<<endl;
    for(int i = 0; i < N + 1; i++)
    {
        cout<< y[i]<<" ";
    }
    cout<<endl<<endl;

    cout<<"The result for method Fourier is: "<<endl;
    for(int i = 0; i < N + 1; i++)
    {
        cout<< prod[i]<<" ";
    }
    cout<<endl<<endl;

    cout<<"The real answer for method Fourier is: "<<endl;
    for(int i = 0; i < N + 1; i++)
    {
        cout<< f[i]<<" ";
    }
    cout<<endl<<endl;
    


  
    return 0;
    
}
