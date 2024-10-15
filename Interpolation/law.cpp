#include "law.hpp"


void law(const int n, const int m, const int a, const int b, vector<double>& nodes)
{

    if(m==1)    // Равноотстоящие
    {
        double h = static_cast<double>(b - a) /(n-1);
        nodes[0] = a;
        for(int i=1;i<n;i++)
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