#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "integration.hpp"
#include <cmath>
using namespace std;
#define pi 3.14159265358979323846



int main(int argc, char **argv)
{   
    int Nx, Ny;
    double Lx, Ly;
    double real = 23. / 45. ;
    printf("Process start \n");
    if(argc>0) cout<<endl;

    Nx = std::stoi(argv[1]);
    Ny = std::stoi(argv[2]);
    Lx = std::stod(argv[3]);
    Ly = std::stod(argv[4]);

    vector<Point> points((Nx + 1) * (Ny + 1));
    vector<Triangle> triangles(2 * Nx * Ny);
    vector<pair<Point, double>> weights(85);
    triangulation(Nx, Ny, Lx, Ly);
    double integral = integrate([](double x, double y) {return x*x*x*x + x*x*y*y + y*y*y*y;}, Nx, points, triangles, weights);
    cout<<"Calculation f(x, y) = x^4 + (xy)^2 + y^4"<<endl;
    cout<<"Numeric Integration = "<< integral << endl;
    cout<<"Real Value          = "<< real     << endl; 


    // int N = Nx;
    // std::ofstream file("result.txt");
    // for (int k = 1; k < 50; ++k) 
    // while(N < 500)
    // {
    //     N += 10;
    //     Nx = N;
    //     Ny = N;
    //     vector<Point> points;
    //     vector<Triangle> triangles;
    //     vector<pair<Point, double>> weights;
    //     triangulation(Nx, Ny, Lx, Ly);
    //     double integral = integrate([](double x, double y) {return x*x*x*x + x*x*y*y + y*y*y*y;}, Nx, points, triangles, weights);
    
        
    //     file
    //     << std::setw(20) 
    //     << std::setprecision(15) 
    //     // << std::log(N) 
    //     << log(N) 

    //     << " "
    //     << std::setw(20) 
    //     << std::setprecision(15) 
    //     // << std::log(std::fabs(integrate([](double x, double y) {return x*x*x*x + x*x*y*y + y*y*y*y;}, points, triangles) - real)) 
    //     << log(std::fabs(integral - real))

    //     << std::endl;
    // }
    // file.close();

    return 0;
    
}
