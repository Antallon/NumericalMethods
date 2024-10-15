#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>                                               //Для setw
#include "solve.hpp"
#include "law.hpp"

#include <math.h>
#include <cmath>
using namespace std;




int main(int argc, char **argv)
{   
    int n,m,a,b;
    ofstream File("result.txt");
    printf("Process start \n");


    n = std::stoi(argv[1]); // string to int
    m = std::stoi(argv[2]);
    a = std::stoi(argv[3]);
    b = std::stoi(argv[4]);

    
    if(argc!=5){
        cerr << "Error: The command line parameters must be 5" << endl;
        return -1;
    }

    if(m<1 || m>3){
        cerr<<"Error: m must be in {1,2,3}"<<endl;
        return -1;
    }
        
    if (!File)
    {
        cerr << "Error: Unable to open file" << endl;
        return -1;
    }

    vector<double> A(n*n);                                       //Вектор для матрицы
    vector<double> y(n);                                         //Значения функции в узлах
    vector<double> y2(n);                                        //Значения функции в узлах (поменяется в процессе solve)
    vector<double> coeff(n);                                     //Вектор коэффициентов полинома
    vector<double> nodes(n);                                     //Вектор узлов

    vector<double> nodes_ans;                                    //Вектор узлов(с добавлением двух равноотстоящих точек между узлами)
    
    law(n,m,a,b,nodes);                                          //Создание узлов по заданному закону

    for(int i = 0;i < n; i++)                                    //Вычисление значений заданной функции в узлах
    {
        y[i] = u(nodes[i]);
        y2[i] = y[i];
    }

    int flag = solve(n,A,y2,coeff,nodes);                         //Решение СЛАУ

    if(flag== -1)
    {
        cerr << "Error: Solving Problem" << endl;
        return -1;
    }

    answer(n, nodes, nodes_ans);                                 //Добавление по две точке между узлами

    File<<std::setw(10)<< " x_k "<<"\t"
    <<std::setw(10)<<" y(x_k) "<<"\t"
    <<std::setw(15)<<" P(x_k) "<<"\t"
    <<std::setw(15)<<" |y(x_k) - P(x_k)| "<<"\t"
    <<std::setw(15)<<" L(x_k) "<<"\t"<<"\t"
    <<std::setw(15)<<" |y(x_k) - L(x_k)| "
    <<endl;


    for(size_t i=0; i < nodes_ans.size(); i++)
    {
        double a = u(nodes_ans[i]),b = P(n, nodes_ans[i], coeff),c = L(n,nodes_ans[i],nodes, y);
        File <<
        std::setw(10) << nodes_ans[i] <<"\t"
        <<std::setw(10)<< a<<"\t"
        <<std::setw(15)<< b << "\t"
        <<std::setw(15)<< std::abs(a - b) << "\t"  
        <<std::setw(15)<< c<<"\t"
        <<std::setw(15)<< std::abs(a - c) << "\t" 
        <<endl;        
    }
    // setw - ширина столбца



    for(size_t i=0; i<coeff.size(); i++)
        cout<<"coef_"<<i<<" = "<< coeff[i] << endl;

    return 0;
    
}