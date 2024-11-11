#include "solve.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;



void computePoints
(
    const int Nx, 
    const int Ny, 
    const double Lx, 
    const double Ly, 
    vector<Point> &points


    
)
{
    double step_x = Lx/Nx;
    double step_y = Ly/Ny;
    int id = 0;
    for(int i = 0; i <=Nx; i++)
    {
        for(int j = 0; j <=Ny; j++)
        {
            points.push_back({id, i * step_x, j * step_y});
            id++;
        }
    }
}

void computeTriangles
(
    const int Nx, 
    const int Ny, 
    vector<Triangle> &triangles
) 
{
    int id = 0;
    for(int i = 0; i < Nx; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            int v1 = i * (Ny + 1) + j;
            int v2 = (i + 1) * (Ny + 1) + j;
            int v3 = i * (Ny + 1) + j + 1;
            int v4 = (i + 1) * (Ny + 1) + j + 1;
            triangles.push_back({id, v1, v2, v4});
            id++;
            triangles.push_back({id, v1, v4, v3});
            id++;
        }
    }
}

void writeTriangulationToFile
(
    const int Nx, 
    const int Ny, 
    const vector<Point> &points, 
    const vector<Triangle> &triangles
) 
{
    ofstream File("triangulation.txt");
    if (!File)
    {
        cerr << "Error: Unable to open file" << endl;
        return;
    }
    File
    << points.size()<<endl
    << triangles.size()<<endl
    << (Nx*(Ny+1) + Ny*(Nx+1) + Nx*Ny) - (2*Nx + 2*Ny)<<endl
    << 2*Nx + 2*Ny
    <<endl;
    for(size_t i = 0; i < points.size(); i++)
    {
        File<<points[i].id<<" "<<points[i].x<<" "<<points[i].y<<endl;
    }
    for(size_t i = 0; i < triangles.size(); i++)
    {
        File<<triangles[i].id<<" "<<triangles[i].p1<<" "<<triangles[i].p2<<" "<<triangles[i].p3<<endl;
    }
}


double integrateTriangle
(
    Point p1, 
    Point p2, 
    Point p3,
    function<double(double, double)> f
)
{
    double f1 = f(p1.x, p1.y);
    double f2 = f(p2.x, p2.y);
    double f3 = f(p3.x, p3.y);
    double area = 0.5 * std::abs((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y));
    return area * (f1 + f2 + f3) / 3.0;
}

double integrate
(
    function<double(double, double)> f,
    const vector<Point> &points,
    const vector<Triangle> &triangles
)
{
    double integral = 0.0;
    for(auto& triangle : triangles)
    {
        integral += integrateTriangle(points[triangle.p1], points[triangle.p2], points[triangle.p3], f);
    }
    return integral;
    
}    
