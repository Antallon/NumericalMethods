#include "solve.hpp"
#include "mnrp.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#define EPS 1e-16
#define eps2 1e-30
using namespace std;





/**
 * @brief Computes a polynomial approximation at a given point.
 * 
 * This function evaluates a polynomial of degree n-1 at a given
 * point x using the provided coefficients. The polynomial is 
 * expressed in the form:
 * 
 * P(x) = coeff[1] * x^0 + coeff[2] * x^1 + ... + coeff[n-1] * x^(n-2)
 * 
 * @param[in] n The number of coefficients, determining the degree of the polynomial.
 * @param[in] x The point at which to evaluate the polynomial.
 * @param[in] coeff A vector containing the polynomial coefficients.
 * 
 * @return The evaluated polynomial value at the point x.
 */
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



/**
 * @brief Constructs a matrix B using polynomial basis functions.
 * 
 * This function populates a matrix B with values computed from the
 * nodes vector using polynomial basis functions. For each node, the
 * function calculates the values for each degree of polynomial up to n-1.
 * Special handling is done for the first two columns to accommodate specific
 * basis functions: the first column alternates between 1.0 and -1.0, and the
 * second column is set to 1.0. For all other columns, the values are calculated
 * as the power of the node value.
 * 
 * @param[in] n The number of nodes and the degree of polynomials plus one.
 * @param[out] B The matrix to store computed polynomial basis values.
 * @param[in] nodes A vector containing the node values used in the computation.
 */
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



/**
 * @brief Iteratively improves the nodes selection for the Valle-Puusen algorithm.
 * 
 * This function iteratively calls VallePuusen to improve the nodes selection for
 * the given approximation. The iteration stops when the VallePuusen algorithm
 * returns 0, indicating that the nodes selection is good enough.
 * 
 * @param[in] n The number of nodes and the degree of polynomials plus one.
 * @param[in] a The left boundary of the interval.
 * @param[in] b The right boundary of the interval.
 * @param[out] B The matrix to store computed polynomial basis values.
 * @param[in,out] nodes A vector containing the node values used in the computation.
 *                      The values are modified by the function.
 * @param[out] coeff Vector containing the coefficients of the best approximation.
 * @param[in] points A vector containing the points at which the function is evaluated.
 */
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




/**
 * @brief Iterative algorithm of Valle-Puusen for finding a grid for the MNRP.
 * 
 * This function implements the iterative algorithm of Valle-Puusen for finding
 * a grid for the MNRP. It takes as input a vector of points, a vector of nodes
 * and a matrix of coefficients. At each iteration, it checks if there is a
 * point in the vector of points that is not in the vector of nodes. If such a
 * point is found, it is added to the vector of nodes and the matrix of
 * coefficients is updated. The algorithm stops when all points in the vector
 * of points are in the vector of nodes.
 * 
 * @param[in] n Number of points in the vector of points.
 * @param[in] B Matrix of coefficients.
 * @param[in] nodes Vector of nodes.
 * @param[in] coeff Vector of coefficients.
 * @param[in] points Vector of points.
 * 
 * @return 0 if the algorithm converged, 1 otherwise.
 */
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







