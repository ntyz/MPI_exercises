#include<iostream>
#include<mpi.h>
#include<ctime>
#include<vector>
#include<cmath>
using namespace std;

#define ull unsigned long long int
#define ll long long int


void get_data(int my_rank, int comm_sz, int *n, double *ndw, double *edw, double *beta)
{
    if (my_rank == 0)
    {
        cout<<"Input n, ndw, edw, beta:"<<endl;
        // cin>>n>>ndw>>edw;
        scanf("%d %lf %lf %lf", n, ndw, edw, beta);
    }
    MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ndw, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(edw, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void get_binary(int *vt, int x, int node_num)
{
    int i = node_num - 1;
    while (x)
    {
        // vt.emplace_back(x & 1);
        vt[i --] = x & 1;
        x >>= 1;
    }
    // while (vt.size() < node_num)  // vt.insert(vt.begin(), 0);
    // vt.emplace_back(0);
    while (i >= 0) vt[i --] = 0;
    
    // return vt;
}

ll get_hamiltonian(double *g, int node_num, int *conf)
{
    ll Hs = 0;
    for (int i = 0; i < node_num; ++ i)
    for (int j = 0; j < node_num; ++ j)
    {
        int confu = 2*int(conf[i]) - 1;
        int confv = 2*int(conf[j]) - 1;
        if (i == j)
        Hs -= g[i * node_num + j] * confu;
        else
        Hs += g[i * node_num + j] * confu * confv;
    }
    return Hs;
}


int main()
{
    int my_rank, comm_sz;
    // clock_t st, ed;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    int n;
    double edw, ndw, beta;
    get_data(my_rank, comm_sz, &n, &ndw, &edw, &beta);    

    // if (my_rank == 0)
    // {
    //     cout<<"Input n, ndw, edw, beta:"<<endl;
    //     cin>>n>>ndw>>edw>>beta;
    // }
    // MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&ndw, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&edw, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int g_sz = n * n;
    double *g = (double *)calloc(g_sz * g_sz, sizeof(double));
    int *conf = (int *)calloc(n, sizeof(int));
    
    for (int i = 0; i < g_sz; ++ i)
    {
        g[i * g_sz + i] = ndw;
        if ((i + 1) % n != 0) g[i * g_sz + i + 1] = g[(i + 1) * g_sz + i] = edw;
        if (i + n < n * n) g[i * g_sz + i + n] = g[(i + n) * g_sz + i] = edw;
    }
    
    for (int i = 0; i < n; ++ i)
    {
        g[i * g_sz + (n - 1) * n + i] = g[((n - 1) * n + i) * g_sz + i] = edw; //列周期

        g[i*n * g_sz + i * n + n - 1] = g[(i * n + n - 1) * g_sz + i * n] = edw; //行周期
    }

    int local_n = pow(2, n) / comm_sz;
    int local_st = my_rank * local_n;
    int local_ed = local_st + local_n - 1;
    double Em = 0.0, psum = 0.0;
    double local_Em = 0.0, local_psum = 0.0;

    double st = MPI_Wtime();
    for (int i = local_st; i <= local_ed; ++ i)
    {
        get_binary(conf, i, n);
        ll Hs = get_hamiltonian(g, n, conf);
        double p = exp(-beta * Hs);
        int mx = 0;
        // for (auto it : conf) mx += 2 * int(it) - 1;
        for (int j = 0; j < n; ++ j) mx += 2 * conf[j] - 1;
        local_psum += p;
        local_Em += p * mx;
    }

    double ed = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&local_Em, &Em, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_psum, &psum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // Em /= psum;

    if (my_rank == 0)
    {
        Em /= psum;
        // cout<<"psum is "<<psum<<endl;
        // cout<<"Em is "<<Em<<endl;
        cout<<"The Ferromagnetism is: "<<Em<<endl;
        cout<<"time = "<<ed - st<<" s"<<endl;

        for (int i = 0; i < g_sz; ++ i)
        {
            for (int j = 0; j < g_sz; ++ j)
            cout<<g[i * g_sz + j]<<" ";
            cout<<endl;
        }
    }
    // else
    // {
    //     cout<<"local_psum is "<<local_psum<<endl;
    // }
    
    free(g); free(conf);
    MPI_Finalize();
    return 0;
}