#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>                                               //Для setw
#include "solve2.hpp"
#include <cmath>

#include <math.h>
#include <cmath>
using namespace std;
#define pi 3.14159265358979323846



int main(int argc,  char **argv)
{   
    if(argc == 2) cout<<endl;
    int N, mIter;
    double p, tau,m, M, q, norma_zero, norma;
    double mult = 1.0;
    ofstream File("result2.txt");
    cout<<"Process start "<<endl<<endl;

    N     = std::stoi(argv[1]);
    mIter = std::stoi(argv[2]);
    p     = std::stod(argv[3]);
   
    if (!File)
    {
        cerr << "Error: Unable to open file" << endl;
        return -1;
    }
    





    vector<double> A((N)*(N));
    vector<double> x(N);
    vector<double> b(N);
    vector<double> prod(N);
    

    m = compute_eigen_value(N, 1, p);
    M = compute_eigen_value(N, N, p);
    tau = 2.0 / (m + M);
    q = (M - m) / (M + m);
    cout<<"q = "<<q<<endl;
    cout<<"m = "<<m<<endl;
    cout<<"M =   "<<M<<endl;
    cout<<"tau = "<< tau<<endl;

    fillmat(N, p, A);
    fillf(N, b);
 
    printMatrix(A,N ,N ,b);
   

    norma_zero = Richardson(x, A, b, tau, N, 1, prod);
    for(int k = 1; k < mIter; k++)
    {
        norma = Richardson(x, A, b, tau, N, k, prod);
        File
        <<setw(20)<< k
        <<setw(20)<< norma
        <<setw(20)<< mult * norma_zero
        << endl;
        mult *= q;
    }
    cout<<"last norma = "<<norma<<endl;
    

    return 0;
    
}
