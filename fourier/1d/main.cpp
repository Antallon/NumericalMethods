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

    try
    {
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
    }

    catch (const std::invalid_argument &e)
    {
        cerr << "Error: Invalid argument. Need enter integers for N" << endl;

        return -1;
    }
    catch(const std::out_of_range&e){
        cerr <<"Error: meaning out of int range"<< endl;
        return -1;
    }
    catch(const std::exception&e){
        cout<<"Error input"<<e.what()<<endl;
        return -1;
    }


    try{

        vector<double> points(N);
        vector<double> p((N-1) * 3 + 1);
        vector<double> coef(N);
        vector<double> diff((N-1)*3 + 1);

        double max_value = -1;


        computePoints(N, points, p);
        coeff(N,points,coef);
        
        for( auto &el : points)
        {
            // cout << "y aprox: " << aprox(N, el, points) << " for x = " << el<< endl;
            // cout << "y real: " << u(el) << " for x = " << el<< endl;
            // cout<< endl;
            File << el <<" " << aprox(N, el, points, coef) << endl;
        }
        File.close();

        for(int i=0;i<diff.size();i++)
        {
            diff[i] = abs(aprox(N, p[i], points, coef) - u(p[i]));
        }

        for(int i=0; i < diff.size(); i++)
        {
            if(diff[i]>max_value)
            {
                max_value = diff[i];
            }

        }
        cout << "Error= " << max_value << endl;

        return 0;
    }
    catch(std::bad_alloc&){
        std::cerr<<"memory error" <<endl;
    }
}
