#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "solve.hpp"
#include <cmath>
using namespace std;
#define pi 3.14159265358979323846



int main(int argc, char **argv)
{   
    int Nx, Ny;
    double Lx, Ly;
    double real = 23. / 45. ;
    ofstream File("result.txt");
    printf("Process start \n");

    Nx = std::stoi(argv[1]);
    Ny = std::stoi(argv[2]);
    Lx = std::stod(argv[3]);
    Ly = std::stod(argv[4]);

    if (!File)
    {
        cerr << "Error: Unable to open file" << endl;
        return -1;
    }
    

    vector<Point> points;
    vector<Triangle> triangles;
    computePoints(Nx, Ny, Lx, Ly, points);
    computeTriangles(Nx, Ny, Lx, Ly, triangles);
    writeTriangulationToFile(Nx, Ny, points, triangles);
    double integral = integrate([](double x, double y) {return x*x*x*x + x*x*y*y + y*y*y*y;}, points, triangles);
    cout<<"Numeric Integration = "<< integral << endl;
    cout<<"Real Value          = "<< real     << endl; 


    int N = Nx;
    std::ofstream file("result.txt");
    for (int k = 1; k < 50; ++k) 
    {
        N += k;
        Nx = N;
        Ny = N;
        vector<Point> points;
        vector<Triangle> triangles;
        computePoints(Nx, Ny, Lx, Ly, points);
        computeTriangles(Nx, Ny, Lx, Ly, triangles);
        file
        << std::setw(20) 
        << std::setprecision(15) 
        << std::log(N) 
        << " "
        << std::setw(20) 
        << std::setprecision(15) 
        << std::log(std::fabs(integrate([](double x, double y) {return x*x*x*x + x*x*y*y + y*y*y*y;}, points, triangles) - real)) 
        << std::endl;
    }
    file.close();

    return 0;
    
}