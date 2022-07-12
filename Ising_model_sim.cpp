// #include<bits/stdc++.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>
// #include "matplotlibcpp.h"
using namespace std;
// namespace plt = matplotlibcpp;
// int main() {
//     plt::plot({1,3,2,4});
//     plt::show();
// }


#define ull unsigned long long int
#define ll long long int
int g[10][10];

void get_chain_g(int node_num, int ndw, int edw) //every weight is same, can be changed to array to set different value for every node & edge
{
    int edi = 0;
    for (int i = 0; i < node_num; ++ i)
    for (int j = 0; j <= i; ++ j)
    {
        if (i == j) g[i][j] = ndw;
        else if (abs(i - j) == 1)
        g[i][j] = g[j][i] = edw;
    }
    g[0][node_num - 1] = g[node_num - 1][0] = edw;
    // return 
}

vector<bool> get_bin(ull x, int node_num)
{
    vector<bool> vt;
    while (x)
    {
        vt.emplace_back(x & 1);
        x >>= 1;
    }
    while (vt.size() < node_num)  // vt.insert(vt.begin(), 0);
    vt.emplace_back(0);
    
    return vt;
}

ll get_hamiltonian(int node_num, vector<bool> conf)
{
    ll Hs = 0;
    for (int i = 0; i < node_num; ++ i)
    for (int j = 0; j < node_num; ++ j)
    {
        int confu = 2*int(conf[i]) - 1;
        int confv = 2*int(conf[j]) - 1;
        if (i == j)
        Hs -= g[i][i] * confu;
        else
        Hs += g[i][j] * confu * confv;
    }
    return Hs;
}

int main()
{
    ull node_num = 3;
    vector<float> betalist;
    for (float beta = 0; beta <= 20; beta += 1) 
        betalist.emplace_back(beta);
    
    get_chain_g(node_num, 1, 1);
    vector<double> ev;
    for (auto beta : betalist)
    {
        double Em = 0.0, psum = 0.0;
        // vector<double> pl;
        for (ull i = 0; i < pow(2UL, node_num); ++ i)
        {
            vector<bool> conf = get_bin(i, node_num);
            ll Hs = get_hamiltonian(node_num, conf);
            double p = exp(-beta * Hs);
            // pl.push_back(p);
            psum += p;
            int mx = 0;
            for (auto it : conf)
                mx += 2 * int(it) - 1;
            // int mx = accumulate(conf.begin(), conf.end(), 0);
            Em += p * mx;
        }
        Em /= psum;
        ev.emplace_back(Em);
    }
    
    for (auto item : ev)
        cout<<item<<" ";
        cout<<endl;
    return 0;
}





// struct node
// {
//     int w
// };

// struct edge
// {
//     int u, v, w;
// };

// void set_node_weight(int node_num, node nd[], int wt)   //传weight数组
// {
//     for (int i = 0; i < node_num; ++ i)
//         nd[i].w = wt;
// }

// void set_edge_weight(int node_num, edge ed[], int wt)
// {
//     for (int i = 0; i < node_num; ++ i)
//         ed[i].w = wt;
// }

// void get_chain_g(int node_num, node nd[], edge ed[])
// {
//     int edi = 0;
//     for (int i = 0; i < node_num; ++ i)
//     {
//         nd[i].w = 1;
//         ed[edi].w = 
//     }
//     for (int i = 0; i < node_num; ++ i)
//     return 
// }

    // int conf[3] = {1, 2, 3};
    // do
    // {
    //     for (int i = 0; i < 3; ++ i)
    //     cout<<conf[i]<<" ";
    //     cout<<endl;
    // } while (next_permutation(conf, conf + 3));