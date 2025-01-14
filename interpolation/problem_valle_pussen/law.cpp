#include "law.hpp"


/**
 * @brief Forms a grid of nodes with given number, type, and bounds.
 * 
 * The type of grid is determined by the parameter m. If m == 1, the nodes are
 * equally spaced. If m == 2, the nodes are the Chebyshev nodes. If m == 3, the
 * nodes are randomly distributed.
 * 
 * @param[in] n Number of nodes
 * @param[in] m Type of grid (1 - equally spaced, 2 - Chebyshev, 3 - random)
 * @param[in] a Lower bound of the grid
 * @param[in] b Upper bound of the grid
 * @param[out] nodes Vector to store the nodes
 */
void law(const int n, const int m, const int a, const int b, vector<double>& nodes)
{

    if(m==1)    // Равноотстоящие
    {
        double h = static_cast<double>(b - a) /(n-1);
        nodes[0] = a;
        for(int i=0;i<n - 1;i++)
        {
            nodes[i] = nodes[i-1] + h;
        }
    }

    if(m==2)   // Узлы Чебышёва
    {
        for(int i=1; i<n+1; i++)
        {
            nodes[i-1] = 0.5*(b-a)*cos(M_PI*(2*i-1)/(2*n)) + 0.5*(a+b);
        }
        reverse(nodes.begin(), nodes.end());
    }

    if(m==3)   // Случайные узлы 
    {
        
    std::mt19937 generator(static_cast<unsigned int>(std::time(0))); //Создание генератора, seed == time(0)

    //time_t time(time_t* timer); time_t - это тип данных. time записывает в *timer текущее время
    //если timer==nullptr, то time просто возвращает текущее время в секундах с начала Unix
    //static_cast явно приводит данные типа time_t к типу unsigned int

    std::uniform_real_distribution<double> distribution(a, b);       //Равномерное распределение

    for (int i = 0; i < n; i++) 
    {
        nodes[i] = distribution(generator);
    }

    sort(nodes.begin(),nodes.end());

    }
}




/**
 * @brief Refines a set of nodes by adding intermediate points between each pair of nodes.
 * 
 * This function takes a vector of nodes and refines it by adding two additional
 * equally spaced points between each pair of consecutive nodes. The refined nodes
 * are stored in the output vector nodes_ans. The spacing between the new points
 * is one third of the interval between the original nodes.
 * 
 * @param[in] n Number of nodes in the input vector
 * @param[in] nodes Vector containing the original nodes
 * @param[out] nodes_ans Vector to store the refined nodes with additional points
 */
void answer(const int n, const vector<double>& nodes, vector<double>& nodes_ans)
{
    for(int i = 0; i < n - 1; i++)
    {
        double h = (nodes[i+1] - nodes[i]) / 3;
        nodes_ans.push_back(nodes[i]);
        nodes_ans.push_back(nodes[i] + h);
        nodes_ans.push_back(nodes[i] + 2*h); 
    }
    nodes_ans.push_back(nodes[n-1]);

}