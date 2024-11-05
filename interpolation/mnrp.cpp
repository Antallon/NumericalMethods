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

    double step = (b - a) / 1000.0;
    coeff[0] = 1.0;

    for(size_t i=0; i < points.size(); i++)
    {
        points[i] = a + i * step;
    }
    // int cou = 0;
    MNRPfull(n, B, nodes);
    vector<double> y(n);
    for(int i = 0;i < n; i++)                                    
    {
        y[i] = u(nodes[i]);
    }
    int flag3 = solve(n, B, y, coeff);
    if(flag3 == -1) cerr << "Error: Solving Problem with Valle-Puusen" << endl;
    // cout<<"h_0 = "<< coeff[0] <<endl;

    int flag = VallePuusen(n, B, nodes, coeff, points);



    // for(int i = 0; i < n; i++) cout<< nodes[i]<< " ";
    // cout<<endl;
    while(true)
    {
        if(flag == 0) break;
        flag = VallePuusen(n, B, nodes, coeff, points);
        // cou++;
    }
    // cout<<"# "<<cou<<endl;
    // cout<<"h_new = "<< coeff[0] <<endl;
    // for(int i = 0; i < n; i++) cout<< nodes[i]<< " ";
    // cout<<endl;

 }




int VallePuusen
(

    const int n, 
    vector<double>& B, 
    vector<double>& nodes, 
    vector<double>& coeff,  
    vector<double>& points
)
{
    int kek = 0;
    vector<double> nodes_temp(n);
    vector<double> B_temp(n*n);
    vector<double> coeff_temp(n);
    vector<double> y_temp(n);

    for(int i = 0; i < n; i++) 
    {
        nodes_temp[i] = nodes[i];
        B_temp[i] = B[i];
        coeff_temp[i] = coeff[i];
    }
    // MNRPfull(n, B, nodes);
    // vector<double> y(n);
    // cout<<"h before = "<< coeff[0]<<endl;


    // for(int i = 0;i < n; i++)                                    
    // {
    //     y[i] = u(nodes[i]);
    // }
    // for(int i = 0; i < n; i++) cout<<nodes[i]<<" ";
    // cout<<endl;


    // int flag3 = solve(n, B, y, coeff);

    // if(flag3 == -1) cerr << "Error: Solving Problem with Valle-Puusen" << endl;
    // cout<<"h after  = "<< coeff[0]<<endl;
    // if(std::abs(coeff[0] - stop) < eps2) return 0;
    // if(std::abs(coeff[0]) < EPS) return 0;
    // cout<<"STOP DIFF = " << stop - coeff[0]<<endl;
    // cout<<"5h = " << 5 * coeff[0]<<endl;



    for( auto&  x : points)
    {   
        for(auto& y : nodes)
        {
            if(std::abs(x - y)< eps2) kek = 1;


        }
        if(kek == 1)
        {
            kek = 0;
            continue;
        }

        // cout<<"x pretendent = "<<x<<endl;
        // cout<<"u(x) = "<<u(x)<<endl;
        // cout<<"mnrp = "<<mnrp(n, x, coeff)<<endl;
        // cout<<"u(x) - mnrp = "<< std::abs(u(x) - mnrp(n, x, coeff))<<endl;
        // cout<<"h = "<<std::abs(  coeff[0])<<endl;
        // cout<<endl; 
//        if(std::abs(u(x) - mnrp(n, x, coeff)) < 5*coeff[0]) continue;
        // if(std::abs(stop - coeff[0])< std::abs(2*coeff[0])) continue;
        if(std::abs(u(x) - mnrp(n, x, coeff)) >  std::abs(coeff[0]))
        {  
            


            // cout<<"/////////////////////////////////////////////////////////////////"<<endl;
            auto it = std::lower_bound(nodes_temp.begin(), nodes_temp.end(), x); // Бинарный поиск,возвращает итератор на первый элемент, который >= x
            // cout<<"it = "<<*it<<endl;
            if(it == nodes_temp.begin())
            {
                if(u(x) * u(nodes_temp[0]) >= 0) nodes_temp[0] = x;
                else
                {
                    nodes_temp.erase(nodes_temp.end() - 1);
                    nodes_temp.insert(nodes_temp.begin(),x);
                }
                // cerr<<"lol1"<<endl;
            }


            else if(it == nodes_temp.end())
            {
                if(u(x) * u(nodes_temp[nodes_temp.size() - 1])>= 0) nodes_temp[nodes_temp.size() - 1] = x;
                else
                {
                    nodes_temp.erase(nodes_temp.begin());
                    nodes_temp.insert(nodes_temp.end(),x);
                }
                // cerr<<"lol2"<<endl;
            }


            else
            {
                double v1 = *(it - 1);
//                double v2 = *it;
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
            // points.erase(points.begin());
            MNRPfull(n, B_temp, nodes_temp);
            for(int i = 0;i < n; i++)                                    
            {
                y_temp[i] = u(nodes_temp[i]);
            }
            int flag3 = solve(n, B_temp, y_temp, coeff_temp);
            if(flag3 == -1) cerr << "Error: Solving Problem with Valle-Puusen" << endl;
            // cout<<"coeff temp = "<<coeff_temp[0]<<endl;
            // cout<<"coeff = "<<coeff[0]<<endl;
            // cout<<"coeff diff = "<<std::abs(coeff_temp[0] - coeff[0])<<endl;
            // cout<<endl;
            if(std::abs(coeff_temp[0] - coeff[0])< 10*std::abs(coeff_temp[0])) continue;
            else
            {
                for(int i = 0; i < n; i++) 
                {
                    nodes[i] = nodes_temp[i];
                    B[i] = B_temp[i];
                    coeff[i] = coeff_temp[i];
                }   
            }
            return 1;
        }

        // points.erase(points.begin());
    }
    return 0;

}







