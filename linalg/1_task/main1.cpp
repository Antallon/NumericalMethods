#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>                                               //Для setw
#include "solve1.hpp"
#include <cmath>

#include <math.h>
#include <cmath>
using namespace std;
#define pi 3.14159265358979323846



int main(int argc,  char **argv)
{   
    if(argc == 2) cout<<endl;
    int N, mIter;
    double p;
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
    





    // vector<double> A((N+1)*(N+1));
    vector<double> A((N)*(N));


    // vector<double> x(N+1);
    vector<long double> y(N+1);
    vector<double> f(N+1);
    vector<double> prod(N);
    


//  Task 1
    fillmat(N, p, A);
    fillf(N, f);
   

    compute_solution_y(N, p, f, y);
    multiply(N, A, y, prod);
    // printMatrix(A, N, N, f);
   


    

















    cout<<"The solution for method Fourier is: "<<endl;
    for(int i = 0; i < N+1 ; i++)
    {
        cout<< y[i]<<" ";
    }
    cout<<endl<<endl;

    cout<<"The result for method Fourier is: "<<endl;
    cout<<"|||||||  ";
    for(int i = 0; i < N ; i++)
    {
        cout<< prod[i]<<" ";
    }
    cout<<endl<<endl;

    cout<<"The real answer for method Fourier is: "<<endl;
    for(int i = 0; i < N+1 ; i++)
    {
        cout<< f[i]<<" ";
    }
    // cout<<endl<<endl;
    // for(int i = 0; i < N ; i++) kek2(A, N, p, i);
    
    return 0;
    
}
