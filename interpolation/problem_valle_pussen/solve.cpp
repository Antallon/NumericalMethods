#include "solve.hpp"
#define EPS 1e-14
#define MAX_ITER 100

double u(double x) {
    //return x * x + x;

    // if (x > 0) {
    //     return x;
    // } 
    // else {
    //     return -x;
    // }

    return std::abs(x);

    // return 1. / (1. + 25 * x * x);
}



double div_diff
(
    const vector<int>& basis, 
    const vector<double>& net, 
    const vector<double>& fun, 
    int start, 
    int end
)
{
    if(end == start)
        return fun[basis[start]];
    return (div_diff(basis, net, fun, start + 1, end) - div_diff(basis, net, fun, start, end - 1)) / (net[basis[end]] - net[basis[start]]);
}



double div_diff_phi
(
    const vector<int>& basis, 
    const vector<double>& net, 
    const vector<double>& phi, 
    int start, 
    int end
)
{
    if(end == start)
        return phi[start];
    return (div_diff_phi(basis, net, phi, start + 1, end) - div_diff_phi(basis, net, phi, start, end - 1)) / (net[basis[end]] - net[basis[start]]);
}


double compute_h
(
    const vector<int>& basis, 
    const vector<double>& net, 
    const vector<double>& fun, 
    const vector<double>& phi, 
    int n
)
{
    return div_diff(basis, net, fun, 0, n + 1)/div_diff_phi(basis, net, phi, 0, n + 1);
}


double compute_polynom_in_point
(
    const vector<double>& fun, 
    const vector<int>&    basis, 
    const vector<double>& net, 
    const vector<double>& phi, 
    double h, 
    double x, 
    int n
)
{
    double result = fun[basis[0]] - h;
    double prod = 1.0;
    double divided_diff_y;
    double divided_diff_phi;

    for (int k = 1; k < n + 1; k++)
    {
        prod *= (x - net[basis[k-1]]);
        divided_diff_y   = div_diff(basis, net, fun, 0, k);
        divided_diff_phi = div_diff_phi(basis, net, phi, 0, k);
        result += (divided_diff_y - h * divided_diff_phi) * prod;
    }
    return result;
}


void compute_polynom_values
(
    vector<double>& poly, 
    const vector<double>& fun, 
    const vector<int>&basis, 
    const vector<double>& nodes, 
    const vector<double>& phi, 
    int n, 
    int N
)
{
    double h = compute_h(basis, nodes, fun, phi, n);
    for(int i = 0; i < N + 1; i++)
    {
        poly[i] = compute_polynom_in_point(fun, basis, nodes, phi, h, nodes[i], n);  
    }     
}


/**
 * @brief Computes the maximum of absolute differences between the function values and the polynomial values on the basis points.
 * @brief Максимальное уклонение на узлах базиса
 *
 * @param basis - vector of indices of the basis points
 * @param fun - vector of function values on the grid
 * @param poly - vector of values of the polynomial on the grid
 * @param n - number of elements in the vectors
 *
 * @return maximum of absolute differences
 */
double compute_h_sigma
(
    const vector<int>&    basis, 
    const vector<double>& fun, 
    const vector<double>& poly, 
    int n
)
{
    double ans = 0.0;
    double temp;
    
    for(int k = 0; k < n + 2; k++)
    {
        temp = abs(fun[basis[k]] - poly[basis[k]]);
        if(temp > ans)
            ans = temp;
    }
    return ans;
}



/**
 * @brief Computes the maximum of absolute differences between elements of two vectors fun and poly.
 * @brief Максимальное уклонение на всей системе узлов
 * 
 * @param fun - vector of function values on the grid
 * @param poly - vector of values of the polynomial on the grid
 * @param N - number of elements in the vectors
 *
 * @return maximum of absolute differences
 */
double compute_phi_sigma
(
    const vector<double>& fun, 
    const vector<double>& poly, 
    int N
)
{
    double ans = 0.0, temp;
    for(int k = 0; k < N + 1; k++)
    {
        temp = abs(fun[k] - poly[k]);
        if(temp > ans)
            ans = temp;
    }

    return ans;
}





void transform_basis
(
    vector<int>& basis, 
    const vector<double>& fun, 
    const vector<double>& poly, 
    double phi_sigma, 
    int n, 
    int N
)
{    
    int num = -1;
    int left = -1;

    for(int i = 0; i < N + 1; i++)
    {
        if(std::abs(std::abs(fun[i] - poly[i]) - phi_sigma) < EPS)
        {
            num = i;
            break;
        } 
    }
   


    if(num < basis[0]) 
    {
        if((fun[num] - poly[num]) * (fun[basis[0]] - poly[basis[0]]) < 0)
        {
            for(int j = n + 1; j > 0; j--)
                basis[j] = basis[j - 1];

        }
        basis[0] = num;
        
    }
    




    else if(num > basis[n + 1]) 
    {
        if((fun[num] - poly[num]) * (fun[basis[n+1]] - poly[basis[n+1]]) < 0)
            for (int j = 0; j < n + 1; j++)
                basis[j] = basis[j + 1];
        basis[n + 1] = num;
        
    }





    else
    {
        for(int i = 0; i < n + 1; i++)
        {
            if((basis[i] < num) && (num < basis[i + 1]))
            {
                left = i;
                break;
            }
        }

        if((fun[num] - poly[num]) * (fun[basis[left]] - poly[basis[left]]) < 0)
            basis[left + 1] = num;
        else
            basis[left] = num;
            return;
    }
    
}


 





 
void vallee_pussen_method
(
    int n, 
    int N,
    vector<int>& basis,
    vector<double>& phi,
    vector<double>& poly,
    const vector<double>& nodes,
    const vector<double>& fun
    
)
{
    int k = 0;
    double h_sigma   = compute_h_sigma(basis, fun, poly, n);
    double phi_sigma = compute_phi_sigma(fun, poly, N);

    // cout << h_sigma << " " << phi_sigma << endl;

    while(k < MAX_ITER && std::abs(h_sigma - phi_sigma)  > (EPS)*1)
    {   
        transform_basis(basis, fun, poly, phi_sigma, n, N);
        compute_polynom_values(poly, fun, basis, nodes, phi, n, N);
        h_sigma = compute_h_sigma(basis, fun, poly, n);
        phi_sigma = compute_phi_sigma(fun, poly, N);
    	// cout << h_sigma << " " << phi_sigma << endl;
        // cout << k << endl;
        k++;
    }
    cout<<setprecision(20)<<"h = "<<h_sigma<<endl;
}
