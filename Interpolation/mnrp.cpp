#include "solve.hpp"
#include "mnrp.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;





double mnrp
(
    const int n, 
    const double x, 
    vector<double>& coeff
)
{
    double sum = 0.0;
    for(int i = 1; i < n; i++)
    {  
        sum += coeff[i] * pow(x, i - 1);
    }
    return sum;
    
}



void MNRPfull
(
    const int n, 
    vector<double>& B, 
    const vector<double>& nodes
)
{
    for(int i=0;i < n; i++)
    {
        for(int j=0; j<n; j++)
        {
            if(j==0) B[i*n + j] = (i % 2 == 0) ? 1.0 : -1.0; // специальный случай для j=0
            else if(j==1) B[i*n + j] = 1.0;
            else
            {
                B[i*n + j] = pow(nodes[i], j-1);
            }
        }
    }
}



void ImproveVallePuusenApproximation
(
    const int n,
    const int a,
    const int b, 
    vector<double>& B, 
    vector<double>& nodes, 
    vector<double>& coeff,  
    vector<double>& points
 )
 {

    double step = (b - a) / 100.0;
    coeff[0] = 1.0;

    for(size_t i=0; i < points.size(); i++)
    {
        points[i] = a + i * step;
    }
    int cou = 0;
    int flag = VallePuusen(n, a, b, B, nodes, coeff, points);
    while(true)
    {
        if(flag == 0) break;
        flag = VallePuusen(n, a, b, B, nodes, coeff, points);
        cou++;
    }
    cout<<"# "<<cou<<endl;

 }




int VallePuusen
(
    const int n, 
    const int a,
    const int b, 
    vector<double>& B, 
    vector<double>& nodes, 
    vector<double>& coeff,  
    vector<double>& points
)
{
    MNRPfull(n, B, nodes);
    vector<double> y(n);
    // cout<<"h before = "<< coeff[0]<<endl;
    double stop = coeff[0];


    for(int i = 0;i < n; i++)                                    
    {
        y[i] = u(nodes[i]);
    }
    // for(int i = 0; i < n; i++) cout<<nodes[i]<<" ";
    // cout<<endl;


    int flag3 = solve(n, B, y, coeff, nodes);

    if(flag3 == -1) cerr << "Error: Solving Problem with Valle-Puusen" << endl;
    // cout<<"h after  = "<< coeff[0]<<endl;
    // if(std::abs(coeff[0] - stop) < eps2) return 0;
    if(std::abs(coeff[0]) < eps2) return 0;
    // cout<<"STOP DIFF = " << stop - coeff[0]<<endl;


    for( auto&  x : points)
    {   
        // cout<<"x pretendent = "<<x<<endl;
        // cout<<"u(x) = "<<u(x)<<endl;
        // cout<<"mnrp = "<<mnrp(n, x, coeff)<<endl;
        // cout<<"u(x) - mnrp = "<< u(x) - mnrp(n, x, coeff)<<endl;
        // cout<<"h = "<<coeff[0]<<endl;
        // cout<<endl; 
        if(std::abs(u(x) - mnrp(n, x, coeff)) < 5*coeff[0]) continue;
        if(std::abs(u(x) - mnrp(n, x, coeff)) > 2.5 * coeff[0])
        {  
            // cout<<"/////////////////////////////////////////////////////////////////"<<endl;
            auto it = std::lower_bound(nodes.begin(), nodes.end(), x); // Бинарный поиск,возвращает итератор на первый элемент, который >= x
            // cout<<"it = "<<*it<<endl;
            if(it == nodes.begin())
            {
                if(u(x) * u(nodes[0]) >= 0) nodes[0] = x;
                else
                {
                    nodes.erase(nodes.end() - 1);
                    nodes.insert(nodes.begin(),x);
                }
                // cerr<<"lol1"<<endl;
            }


            else if(it == nodes.end())
            {
                if(u(x) * u(nodes[nodes.size() - 1])>= 0) nodes[nodes.size() - 1] = x;
                else
                {
                    nodes.erase(nodes.begin());
                    nodes.insert(nodes.end(),x);
                }
                // cerr<<"lol2"<<endl;
            }


            else
            {
                double v1 = *(it - 1);
                double v2 = *it;
                // cout<<"v1 = "<<v1<<endl;
                // cout<<"v2 = "<<v2<<endl;
                // for(int i = 0; i < n; i++) cout<<nodes[i]<<" ";
                // cout<<endl;
                if(u(x) * u(v1) >= 0) *(it - 1) = x;
                else *it = x;
                // cout<<"v1 = "<<v1<<endl;
                // cout<<"v2 = "<<v2<<endl;
                // for(int i = 0; i < n; i++) cout<<nodes[i]<<" ";
                // cout<<endl;

                // cerr<<"lol3"<<endl;
            }
            points.erase(points.begin());
            return 1;
        }
        points.erase(points.begin());
    }
    return 0;

}







