#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>                                               //Для setw
#include "solve.hpp"
#include "law.hpp"
#include "polynom.hpp"
#include "mnrp.hpp"

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
    vector<double> B(n*n);
    vector<double> y(n);                                         //Значения функции в узлах
    vector<double> y2(n);                                        //Значения функции в узлах (поменяется в процессе solve)
    vector<double> y3(n);
    vector<double> coeff(n);                                     //Вектор коэффициентов полинома
    vector<double> coeff_mnrp(n);
    vector<double> coeff_vp(n);
    vector<double> nodes(n);                                     //Вектор узлов
    vector<double> nodes_vp(n);
    vector<double> points(1000);                                 //Вектор точек для Валле-Пуссена

    vector<double> nodes_ans;                                    //Вектор узлов(с добавлением двух равноотстоящих точек между узлами)
    
    law(n,m,a,b,nodes);                                          //Создание узлов по заданному закону

    for(int i = 0;i < n; i++)                                    //Вычисление значений заданной функции в узлах
    {
        y[i] = u(nodes[i]);
        y2[i] = y[i];
        y3[i] = y[i];
        nodes_vp[i] = nodes[i];
    }

    Pmatfull(n, A, nodes);                                       //Создание матрицы для полинома
    MNRPfull(n, B, nodes);                                       //Создание матрицы для МНРП
    int flag1 = solve(n,A,y2,coeff);                             //Решение СЛАУ для полинома
    int flag2 = solve(n,B,y3,coeff_mnrp);                        //Решение СЛАУ для МНРП
    ImproveVallePuusenApproximation(n, a, b, B, nodes_vp, coeff_vp, points);
    if(flag1 == -1)
    {
        cerr << "Error: Solving Problem with Polynomial" << endl;
        return -1;
    }
    if(flag2 == -1)
    {
        cerr << "Error: Solving Problem with MNRP" << endl;
        return -1;
    }

    answer(n, nodes, nodes_ans);                                 //Добавление по две точке между узлами

    File
    <<std::setw(10)<<" x_k "                 <<"\t"
    <<std::setw(10)<<" y(x_k) "              <<"\t"
    <<std::setw(15)<<" P(x_k) "              <<"\t"
    <<std::setw(15)<<" |y(x_k) - P(x_k)| "   <<"\t"
    <<std::setw(10)<<" L(x_k) "              <<"\t"
    <<std::setw(15)<<" |y(x_k) - L(x_k)| "   <<"\t"
    <<std::setw(10)<<" MNRP(x_k) "    <<"\t"
    <<std::setw(15)<<" |y(x_k) - MNRP(x_k)| "<<"\t"
    <<std::setw(10)<<" VP(x_k) "      <<"\t"
    <<std::setw(15)<<" |y(x_k) - VP(x_k)|   "<<"\t"
    <<endl;


    for(size_t i=0; i < nodes_ans.size(); i++)
    {
        double a = u(nodes_ans[i]),b = P(n, nodes_ans[i], coeff),c = L(n,nodes_ans[i],nodes, y);
        File
        <<std::setw(10)<< nodes_ans[i]                                   <<"\t"
        <<std::setw(10)<< a                                              <<"\t"
        <<std::setw(15)<< b                                              <<"\t"
        <<std::setw(15)<< std::abs(a - b)                                <<"\t"  
        <<std::setw(15)<< c                                              <<"\t"
        <<std::setw(15)<< std::abs(a - c)                                <<"\t" 
        <<std::setw(15)<< mnrp(n, nodes_ans[i], coeff_mnrp)              <<"\t"
        <<std::setw(15)<< std::abs(a - mnrp(n, nodes_ans[i], coeff_mnrp))<<"\t"
        <<std::setw(15)<< mnrp(n, nodes_ans[i], coeff_vp)                <<"\t"
        <<std::setw(15)<< std::abs(a - mnrp(n, nodes_ans[i], coeff_vp))  <<"\t"
        <<endl;        
    }
    // setw - ширина столбца



    // for(size_t i=0; i<coeff.size(); i++)
    //     cout<<"coef_"<<i<<" = "<< coeff[i] << endl;
//    cout<<"h = "<<coeff_mnrp[0]<<endl;
    return 0;
    
}
