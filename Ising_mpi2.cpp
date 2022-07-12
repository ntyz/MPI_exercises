#include<iostream>
#include<mpi.h>
#include<ctime>
#include<vector>
#include<cmath>
using namespace std;

#define ull unsigned long long int
#define ll long long int

void get_data(int my_rank, int comm_sz, int *n, int *m, double *beta)
{
    if (my_rank == 0)
    {
        cout<<"Input n, m, beta:"<<endl;
        // std::cin>>n>>m;  
        scanf("%d %d %lf", n, m, beta);
    }
    MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

ll get_hamiltonian(int node_id, int n, int m, int *conf)
{
    ll Hs = 0;
    int rid = node_id / m, cid = node_id % m;
    int dwr = (rid + 1) % n, dwc = cid;
    int upr = (rid - 1 + n) % n, upc = cid;
    int rtr = rid, rtc = (cid + 1) % m;
    int ler = rid, lec = (cid - 1 + m) % m;

    int uid = upr * m + upc;
    int did = dwr * m + dwc;
    int leid = ler * m + lec; 
    int rtid = rtr * m + rtc;
    Hs = conf[node_id] * conf[uid];
    if (uid != did)
        Hs += conf[node_id] * conf[did];
    
    Hs += conf[node_id] * conf[rtid];
    if (leid != rtid)
        Hs += conf[node_id] * conf[leid];
    // cout<<"node id is: "<<node_id<<" Hs is: "<<Hs<<endl;
    return Hs;
}

int main()
{
    int wd_rank, wd_sz;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &wd_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &wd_rank);

    int n, m;
    double beta;
    get_data(wd_rank, wd_sz, &n, &m, &beta);   
    int *conf = (int *)calloc(n * m, sizeof(int));
    for (int i = 0; i < n * m; ++ i)
    conf[i] = (i & 1) == 0 ? 1 : -1;

    // int color = wd_rank / m;
    // MPI_Comm sub_comm;
    // int sub_rank, sub_size;
    // MPI_Comm_split(MPI_COMM_WORLD, color, wd_rank, &sub_comm);
    // MPI_Comm_rank(sub_comm, &sub_rank);
    // MPI_Comm_size(sub_comm, &sub_size);

    double Hs = 0.0;
    double st = MPI_Wtime();
    double local_Hs = get_hamiltonian(wd_rank, n, m, conf);
    // cout<<"rank is: "<<wd_rank<<" local_Hs is: "<<local_Hs<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&local_Hs, &Hs, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double ed = MPI_Wtime();

    if (wd_rank == 0)
    {
        int mx = 0;
        for (int j = 0; j < n * m; ++ j) mx += 2 * conf[j] - 1;
        Hs /= 2;
        cout<<"Hs is: "<<Hs<<endl;
        double local_p = exp(-beta * Hs);
        double local_Em = local_p * mx;
        cout<<"The local Ferromagnetism is: "<<local_Em<<endl;
        cout<<"time = "<<ed - st<<" s"<<endl;
    }
    // else
    // {
    //     cout<<"local_psum is "<<local_psum<<endl;
    // }
    
    free(conf);

    MPI_Finalize();
    return 0;
}