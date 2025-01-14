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
