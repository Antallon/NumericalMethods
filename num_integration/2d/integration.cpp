#include "integration.hpp"

#include "integration.hpp"
void triangulation(int Nx, int Ny, double Lx, double Ly)
{
    double hx = Lx / Nx;
    double hy = Ly / Ny;
    int numNodes = (Nx + 1) * (Ny + 1);
    int numTriangles = 2 * Nx * Ny; 
    int numInnerEdges = (3 * Nx + 2) * Ny;
    int numBoundaryEdges = 2 * (Nx + Ny);
    ofstream out("triangulation.txt");

    out << numNodes << endl;
    out << numTriangles << endl;
    out << numInnerEdges << endl;
    out << numBoundaryEdges << endl;

    for(int i = 0; i < Nx + 1; i++)
        for(int j = 0; j < Ny + 1; j++) 
            out << j + i * (Ny + 1) + 1 << ": " << fixed << setprecision(15) << i * hx << " " << j * hy << endl;


    for(int i = 0; i < Nx; i++)
        for(int j = 0; j < Ny; j++) 
        {
            out << 2 * (j + i * Ny) + 1 << ": " << j + i * (Ny + 1) + 1 << " " << j + i * (Ny + 1) + 2 << " " << j + (i + 1) * (Ny + 1) + 1 << endl;
            out << 2 * (j + i * Ny) + 2 << ": " << j + (i + 1) * (Ny + 1) + 1 << " " << j + (i + 1) * (Ny + 1) + 2 << " " << j + i * (Ny + 1) + 2 << endl;
        }

  
    int k = 1;

    for(int i = 0; i < Nx; i++)
        for(int j = 0; j < Ny; j++)
            out << k++ << ": " << j + (i + 1) * (Ny + 1) + 1 << " " << j + i * (Ny + 1) + 2 << endl;

    for(int i = 1; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            out << k++ << ": " << j + i * (Ny + 1) + 1 << " " << j + (i + 1) * (Ny + 1) + 1 << endl;

    for(int i = 0; i < Nx; i++)
        for (int j = 1; j < Ny; j++)
            out << k++ << ": " << j + i * (Ny + 1) + 1 << " " << j + i * (Ny + 1) + 2 << endl;

    k = 1;  
    for (int i = 0; i < Nx; i++)
        out << k++ << ": " << i + 1 << " " << i + 2 << endl;

    for (int j = 0; j < Ny; j++)
        out << k++ << ": " << (j + 1) * (Nx + 1) << " " << (j + 2) * (Nx + 1) << endl;

    for (int i = Nx; i > 0; i--)
        out << k++ << ": " << Ny * (Nx + 1) + i + 1 << " " << Ny * (Nx + 1) + i + 2 << endl;
    
    for (int j = Ny; j > 0; j--)
        out << k++ << ": " << j * (Nx + 1) + 1 << " " << (j + 1) * (Nx + 1) + 1 << endl;

    out.close();
}




void get_vertices_and_triangles_from_file
(
    vector<Point>& vertices,
    vector<Triangle>& triangles
)
{
    int index;
    char two_points;
    ifstream file("triangulation.txt");
    int num_vert, num_trian, num_inter, num_bound;
    file   >>
    num_vert  >>
    num_trian >>
    num_inter >>
    num_bound;


    for(size_t i = 0; i < vertices.size(); i++)
    {
        
        file >> index >> two_points >> vertices[i].x >> vertices[i].y;
        
    }

    for(size_t i = 0; i < triangles.size(); i++) 
    {
        int index;
        char two_points;
        file >> index >> two_points >> triangles[i].v1 >> triangles[i].v2 >> triangles[i].v3;

        auto &v1 = vertices[triangles[i].v1 - 1];
        auto &v2 = vertices[triangles[i].v2 - 1];
        auto &v3 = vertices[triangles[i].v3 - 1];

        triangles[i].p1 = {v1.x, v1.y};
        triangles[i].p2 = {v2.x, v2.y};
        triangles[i].p3 = {v3.x, v3.y};
    }
    file.close();
}



void get_weights_from_file
(
    vector<pair<Point, double>>& weights
)
{
    ifstream file("weights.txt");

    double x, y, weight;
    while(file >> x >> y >> weight)
    {
        Point point = {x, y};
        weights.push_back({point, weight});
    }

    file.close();
}




double integrate
(
    function<double(double, double)> f,
    int N,
    vector<Point> &vertices,
    vector<Triangle> &triangles,
    vector<pair<Point, double>> &weights
)
{
    int num_vert, num_trian, num_inter, num_bound;
    double integral = 0.0;
    Point A, B, C;
    ifstream filein("triangulation.txt");

    
    get_vertices_and_triangles_from_file(vertices, triangles);
    get_weights_from_file (weights);


    for(size_t i = 0; i < triangles.size(); i++)
    {
        A = triangles[i].p1;
        B = triangles[i].p2;
        C = triangles[i].p3;
        for(size_t j = 0; j < weights.size(); j++)
        {
            Point p;
            p.x = A.x + weights[j].first.x * (B.x - A.x) + weights[j].first.y * (C.x - A.x);
            p.y = A.y + weights[j].first.x * (B.y - A.y) + weights[j].first.y * (C.y - A.y);
            double area = 0.5 * abs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
            // double area = 0.5 * std::abs((B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y));

            integral += 2 * area * f(p.x, p.y) * weights[j].second;
        }
        // double area = 0.5 * abs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
        // double area = 0.5 * abs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));

        // integral += area * (f(A.x, A.y) + f(B.x, B.y) + f(C.x, C.y)) / 3.0;
    }
    filein.close();
    return integral;
}



