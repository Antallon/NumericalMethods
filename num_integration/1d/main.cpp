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
    if(argc == 1) cout<<endl;
    int N;
    ofstream File("result.txt");
    printf("Process start \n");

    N = std::stoi(argv[1]);
    // a = std::stoi(argv[3]);
    // b = std::stoi(argv[4]);

    // if(argc!=5){
    //     cerr << "Error: The command line parameters must be 5" << endl;
    //     return -1;
    // }  


    if (!File)
    {
        cerr << "Error: Unable to open file" << endl;
        return -1;
    }
    cout<<endl<<endl;
    cout<<"f(x) = 1, x in [1.0, 1.1]:"<<endl;
    cout<<"Simpson Integral = "<<Integral_Simpson(1.0, 1.1, [](double x){return x*0 + 1;})<<endl;
    cout<<"Gauss Integral   = "<<Integral_Gauss  (1.0, 1.1, [](double x){return x*0 + 1;})<<endl;
    cout<<"Real  Answer     = 0.1"<<endl;
    cout<<endl<<endl;

    cout<<"f(x) = x, x in [1.0, 1.1]:"<<endl;
    cout<<"Simpson Integral = "<<Integral_Simpson(1.0, 1.1, [](double x){return x;})<<endl;
    cout<<"Gauss Integral   = "<<Integral_Gauss  (1.0, 1.1, [](double x){return x;})<<endl;
    cout<<"Real  Answer     = 0.105"<<endl;
    cout<<endl<<endl;

    cout<<"f(x) = x^2, x in [1.0, 1.1]:"<<endl;
    cout<<"Simpson Integral = "<<Integral_Simpson(1.0, 1.1, [](double x){return x*x;})<<endl;
    cout<<"Gauss Integral   = "<<Integral_Gauss  (1.0, 1.1, [](double x){return x*x;})<<endl;
    cout<<"Real  Answer     = 0.110333"<<endl;
    cout<<endl<<endl;

    cout<<"f(x) = x^3, x in [1.0, 1.1]:"<<endl;
    cout<<"Simpson Integral = "<<Integral_Simpson(1.0, 1.1, [](double x){return x*x*x;})<<endl;
    cout<<"Gauss Integral   = "<<Integral_Gauss  (1.0, 1.1, [](double x){return x*x*x;})<<endl;
    cout<<"Real  Answer     = 0.116025"<<endl;
    cout<<endl<<endl;

    cout<<"f(x) = x^5, x in [1.0, 1.1]:"<<endl;
    cout<<"Simpson Integral = "<<Integral_Simpson(1.0, 1.1, [](double x){return x*x*x*x*x;})<<endl;
    cout<<"Gauss Integral   = "<<Integral_Gauss  (1.0, 1.1, [](double x){return x*x*x*x*x;})<<endl;
    cout<<"Real  Answer     = 0.128594"<<endl;
    cout<<endl<<endl;

    cout<<"f(x) = x^9, x in [1.0, 1.1]:"<<endl;
    cout<<"Simpson Integral = "<<Integral_Simpson(1.0, 1.1, [](double x){return x*x*x*x*x*x*x*x*x;})<<endl;
    cout<<"Gauss Integral   = "<<Integral_Gauss  (1.0, 1.1, [](double x){return x*x*x*x*x*x*x*x*x;})<<endl;
    cout<<"Real  Answer     = 0.159374"<<endl;
    cout<<endl<<endl;

    cout<<"///////////////////////////////////////////////////////////////////"<<endl;

    cout<<endl<<endl;
    cout<<"f(x) = cos(100x), x in [0,pi]:"<<endl;
    cout<<"Simpson Integral Composite = "<<Integral_Simpson_Composite(N, 0.0, pi, [](double x){return cos(100.0*x);})<<endl;
    cout<<"Simpson Integral Composite = "<<Integral_Gauss_Composite  (N, 0.0, pi, [](double x){return cos(100.0*x);})<<endl;
    cout<<"Real  Answer               = 0"<<endl;
    cout<<endl<<endl;

//    cout<<"f(x) = exp(-1000x), x in [0,1]:"<<endl;
//    cout<<"Simpson Integral Composite = "<<Integral_Simpson_Composite(N, 0.0, 1.0, [](double x){return exp(-1000.0*x);})<<endl;
//    cout<<"Simpson Integral Composite = "<<Integral_Gauss_Composite  (N, 0.0, 1.0, [](double x){return exp(-1000.0*x);})<<endl;
//    cout<<"Real  Answer               ~ 10^(-3)"<<endl;
//    cout<<endl<<endl;

    cout<<"f(x) = (1 - x^2)^(-0.5), x in [-1,1]:"<<endl;
    cout<<"Simpson Integral Composite = "<<Integral_Simpson_Composite(N, -1.0, 1.0, [](double x){return 1.0 / sqrt(1.0 - x * x + 1e-6);})<<endl;
    cout<<"Simpson Integral Composite = "<<Integral_Gauss_Composite  (N, -1.0, 1.0, [](double x){return 1.0 / sqrt(1.0 - x * x);})<<endl;
    cout<<"Real  Answer               = pi"<<endl;
    cout<<endl<<endl;

    cout<<"f(x) = e^x, x in [0,20]:"<<endl;
    cout<<"Simpson Integral Composite = "<<Integral_Simpson_Composite(N, 0.0, 20.0, [](double x){return exp(x);})<<endl;
    cout<<"Simpson Integral Composite = "<<Integral_Gauss_Composite  (N, 0.0, 20.0, [](double x){return exp(x);})<<endl;
    cout<<"Real  Answer               = 485 165 194.41"<<endl;
    cout<<endl<<endl;

    ofstream file("result.txt");
    while(N < 100000) {

        
    
        N += 1000;

        double LS1 = log(fabs(Integral_Simpson_Composite(N, 0.0001, 1.0, [](double x){return 1/sqrt(x);}) - 2*(sqrt(1.0) - sqrt(0.0001))));
        double LG1 = log(fabs(Integral_Gauss_Composite  (N, 0.0001, 1.0, [](double x){return 1/sqrt(x);}) - 2*(sqrt(1.0) - sqrt(0.0001))));
        double LS2 = log(fabs(Integral_Simpson_Composite(N, 0.001,  1.0, [](double x){return 1/sqrt(x);}) - 2*(sqrt(1.0) - sqrt(0.001))));
        double LG2 = log(fabs(Integral_Gauss_Composite  (N, 0.001,  1.0, [](double x){return 1/sqrt(x);}) - 2*(sqrt(1.0) - sqrt(0.001))));
        double LS3 = log(fabs(Integral_Simpson_Composite(N, 0.00001,   1.0, [](double x){return 1/sqrt(x);}) - 2*(sqrt(1.0) - sqrt(0.00001))));
        double LG3 = log(fabs(Integral_Gauss_Composite  (N, 0.00001,   1.0, [](double x){return 1/sqrt(x);}) - 2*(sqrt(1.0) - sqrt(0.00001))));
        file 
        << setw(15) << log(N)<< " "
        << setw(15)  << LS1 << " "
        << setw(15)  << LG1 << " "
        << setw(15)  << LS2 << " "
        << setw(15)  << LG2 << " "
        << setw(15)  << LS3 << " "
        << setw(15)  << LG3 << " "
        <<endl;
    }


    


  


   
  
    return 0;
    
}
