#include <iostream>
#include <fstream>
#include <vector>
#include "solve.hpp"

#include <math.h>
#include <cmath>
using namespace std;


int main(int argc, char **argv)
{   
    int N;
    ofstream File("result.txt");
    printf("Process start \n");


    N = std::stoi(argv[1]); // string to int
    if(N < 2){
        cerr << "Error: N must be >= 2" << endl;
        return -1;
    }
    if(argc!=2){
        cerr << "Error: The command line parameters must be 2" << endl;
        return -1;
    }
        
    if (!File)
    {
        cerr << "Error: Unable to open file" << endl;
        return -1;
    }
    

    


    
    vector<double> U(N * N);
    vector<double> D(N * N);
    vector<double> C(N * N);

    vector<double> points(N);
    vector<double> p((N-1) * 3 + 1);
    vector<double> coef(N);
    vector<double> diff;

    double max_value = -1;
    double lol1 = 0;

    computePoints(N, points, p);
    Ufull(N, U, points);
    Cfull(N, points, C, U);
    Dfull(N, points, C, D);
        
    for( auto &elx : points)
    {
        for(auto &ely : points)
        {
            File << elx <<" " << ely <<" "<< aprox2(N, elx, ely, points, U, D) << endl;
        }
    }
    File.close();

    for(auto &elx : p)
    {
        for(auto &ely:p)
        {
            diff.push_back(abs(aprox2(N,elx,ely,points,U,D) - u(elx,ely)));
        }
    }

    for(int i=0; i < diff.size(); i++)
    {
        if(diff[i]>max_value)
        {
            max_value = diff[i];
        }

    }
        cout << "Error = " << max_value << endl;

     

    return 0;
    
}
