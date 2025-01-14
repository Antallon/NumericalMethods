#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>                                               //Для setw
#include "solve3.hpp"
#include <cmath>

#include <math.h>
#include <cmath>
using namespace std;
#define pi 3.14159265358979323846



int main(int argc,  char **argv)
{   
    if(argc == 2) cout<<endl;
    int N, mIter;
    double p, q, norma_zero;
    double tau = 1.0;
    // double tau = 0.1;
    double mult = 1.0;
    ofstream File("result3.txt");
    cout<<"Process start "<<endl<<endl;

    N     = std::stoi(argv[1]);
    mIter = std::stoi(argv[2]);
    p     = std::stod(argv[3]);
   
    if (!File)
    {
        cerr << "Error: Unable to open file" << endl;
        return -1;
    }
     





    vector<double> A(N*N);
    vector<double> B(N*N);
    vector<double> x(N+1);       //x_1,....,x_N
    vector<double> b(N+1);       //b_1,....,b_N
    vector<double> prod(N+1);    //(Ax)_1,....,(Ax)_N     
    vector<double> prod1(N+1); // y_0,....,y_N
    vector<double> prod2(N+1); // y_0,....,y_N
    // vector<double> x_random(N+1);
    // x_random[0] = 0.0;
    // for(int i = 1; i < N+1; i++)
    //     x_random[i] = cos(M_PI*1.0/double(i));


    fillMat(N, p, A, true);
    fillMat(N, p, B, false);

    fillf(N, b);
    // fillf2(B, N, b, x_random);
   



   for(int i = 0; i < N+1; i++) x[i] = 0.0;

    // printMatrix(A, N, N, b);
    // printMatrix(B, N, N, b);
    q = Findq(A, N);
    cout<<"q = "<< q<<endl;
    norma_zero = BSolver(x, A, B, b, tau, N, 0, prod, prod1,prod2, p);

    for(int i = 0; i < N+1; i++) x[i] = 0.0;

    for (int k = 0; k < mIter; k++)
    {
        double norm = BSolver(x, A, B, b, tau, N, k, prod, prod1,prod2, p);
        File
        <<setw(20)<< k
        <<setw(20)<< norm
        <<setw(20)<< mult * norma_zero
        << endl;
        mult *= q;
    }
    
    return 0;
    
}
