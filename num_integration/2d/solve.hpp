#pragma once
#include <cmath>
#include <vector>
#include <functional>

using namespace std;

struct Point
{
    int id;
    double x, y;
    
};
struct Triangle
{
    int id;
    int p1, p2, p3;
};

void computePoints(const int Nx,const int Ny, const double Lx, const double Ly, vector<Point> &points);
void computeTriangles(const int Nx, const int Ny, const double Lx, const double Ly, vector<Triangle> &triangles);
void writeTriangulationToFile(const int Nx, const int Ny, const vector<Point> &points, const vector<Triangle> &triangles);
double integrateTriangle(Point p1, Point p2, Point p3,function<double(double, double)> f);
double integrate(function<double(double, double)> f,const vector<Point> &points,const vector<Triangle> &triangles);