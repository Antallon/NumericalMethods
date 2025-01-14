#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>                                               //Для setw
#include "solve.hpp"
#include "law.hpp"
#include <cmath>
using namespace std;




int main(int argc, char **argv)
{   
    int n,m,a,b;
    int N = 1000;
    
    double ans = 0.0;
    ofstream File("result.txt");
    printf("Process start \n");


    n = std::stoi(argv[1]);
    m = std::stoi(argv[2]);
    a = std::stoi(argv[3]);
    b = std::stoi(argv[4]);
    int h = N / (n + 1);
    


    vector<double> poly(N + 1);                       //Значения полинома на решетке
    vector<int>    basis(n + 2);                      //вектор индексов - наш базис
    vector<double> fun(N + 1);                        //Значения функции на решетке                  
    vector<double> nodes(N + 1);                                          
    vector<double> phi(n + 2);
    
    law(N + 1,m,a,b,nodes);                                          //Создание узлов по заданному закону
    
    for(int i = 0; i < N + 1; i++)                                    
    {
        fun[i] = u(nodes[i]);
    }

    
    for(int k = 0; k < n + 2; k++)
    {
        basis[k] = k*h;
        phi[k] = pow(-1, k);
    }

    compute_polynom_values(poly, fun, basis, nodes, phi, n, N);       //Заполнение вектора poly значениями полинома на решетке
    vallee_pussen_method(n, N, basis, phi, poly, nodes, fun);         //Вычисление значения функции на решетке


    for(int i = 0; i < N + 1; i++)
    {
        ans = std::max(ans, std::abs(poly[i] - fun[i]));
    }

    for(size_t i = 0; i < fun.size(); i++){

        File 
        << std::setprecision(15)
        << nodes[i] 
        << " " 
        << fun[i] 
        << " "
        << poly[i] 
        << " " 
        << std::abs(fun[i] -  poly[i]) 
        << endl;
    }
    cout<<"kek"<<endl;
    return 0;
    
}
