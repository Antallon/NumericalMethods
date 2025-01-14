#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <functional>

using namespace std;
#define EPS 1e-5

struct Point 
{
    double x, y;
};

struct Triangle 
{
    int v1, v2, v3;
    Point p1, p2, p3;
};

double integrate(function<double(double, double)> f,int N,vector<Point> &vertices,vector<Triangle> &triangles,vector<pair<Point, double>> &w);
void get_vertices_and_triangles_from_file(vector<Point>& vertices,vector<Triangle>& triangles);
// void get_triangle_from_file(const vector<Point> &vertices, vector<Triangle>& triangles);
void get_weights_from_file(vector<pair<Point, double>>& weights);
void triangulation(int Nx, int Ny, double Lx, double Ly);

