#include<iostream>
#include<mpi.h>
#include<ctime>
#include<vector>
#include<cmath>
#include<unistd.h>
using namespace std;

#define ull unsigned long long int
#define ll long long int

void get_data(int my_rank, int comm_sz, int *r, int *c, double *beta)
{
    if (my_rank == 0)
    {
        cout<<"Input r, c, beta:"<<endl;
        // std::cin>>n>>m;
        scanf("%d %d %lf", r, c, beta);
    }
    MPI_Bcast(r, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(c, 1, MPI_INT, 0, MPI_COMM_WORLD);
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

double get_hamiltonian(int sub_rank, int n, int m, int *conf)
{
    double Hs = 0;
    for (int i = 0; i < m; ++ i)
    {
        int node_id = sub_rank * m + i;
        int rid = node_id / m, cid = node_id % m;

        int dwr = (rid + 1) % n, dwc = cid;
        int upr = (rid - 1 + n) % n, upc = cid;
        int rtr = rid, rtc = (cid + 1) % m;
        int ler = rid, lec = (cid - 1 + m) % m;

        int uid = upr * m + upc;
        int did = dwr * m + dwc;
        int leid = ler * m + lec; 
        int rtid = rtr * m + rtc;

        int confu = 2*int(conf[uid]) - 1;
        int confd = 2*int(conf[did]) - 1;
        int confl = 2*int(conf[leid]) - 1;
        int confr = 2*int(conf[rtid]) - 1;
        int confi = 2*int(conf[node_id]) - 1;

        if (node_id != uid)
            Hs +=  confi * confu;
        if (uid != did && node_id != did)
            Hs += confi * confd;
        
        if (node_id != rtid)
            Hs += confi * confr;
        if (leid != rtid && node_id != leid)
            Hs += confi * confl;
        // cout<<"node id is: "<<node_id<<" Hs is: "<<Hs<<endl;
    }
    Hs /= 2.0;
    // cout<<"sub_rank is: "<<sub_rank<<" Hs is: "<<Hs<<endl;
    return Hs;
}


int main()
{

    // int i = 0;
    // while (i == 0)
    // {
    //     sleep(5);
    // }

    int wd_rank, wd_sz;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &wd_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &wd_rank);

    int r = 7, c = 7;
    double beta = 1;
    // get_data(wd_rank, wd_sz, &r, &c, &beta);   
    int *conf = (int *)calloc(r * c, sizeof(int));
    // double *beta = (double *)calloc(beta_len, sizeof(int));

    ll node_num = r * c;
    ll conf_num = pow(2, node_num);
    if (wd_sz % r != 0 || conf_num % (wd_sz / r) != 0)
    {
        if (wd_rank == 0)
        {
            cout<<"the value wd_sz % r should be 0"<<endl;
            cout<<"the value conf_num % (wd_sz / r) should be 0"<<endl;
        }
        exit (0);
    }

    clock_t st, ed;

    int group_num = wd_sz / r;
    int local_n = conf_num / group_num;
    // printf("group num is %d , local_n is %d\n", group_num, local_n);

    int color = wd_rank / (wd_sz / group_num); 
    MPI_Comm sub_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wd_rank, &sub_comm);
    int sub_rank, sub_size;
    MPI_Comm_rank(sub_comm, &sub_rank);
    MPI_Comm_size(sub_comm, &sub_size);

    int local_gst = wd_rank / r * local_n;
    int local_ged = local_gst + local_n - 1;
    double Em = 0.0;    //, psum = 0.0;
    double local_Em = 0.0;   //, local_psum = 0.0;
    double Hs = 0, mx = 0;

    // printf("WORLD_RANK is: %d, local_gst is: %d\n", wd_rank, local_gst);
    st = clock(); 
    for (int i = local_gst; i <= local_ged; ++ i)
    {
        Hs = 0; mx = 0;
        get_binary(conf, i, node_num);
        double local_Hs = get_hamiltonian(sub_rank, r, c, conf);
        MPI_Barrier(sub_comm);
        MPI_Reduce(&local_Hs, &Hs, 1, MPI_DOUBLE, MPI_SUM, 0, sub_comm);

        if (sub_rank == 0)
        {
            double p = exp(-beta * Hs);
            for (int j = 0; j < node_num; ++ j) mx += 2 * conf[j] - 1;
            mx /= node_num;
            local_Em += mx * p;
        } 

        // MPI_Barrier(MPI_COMM_WORLD);
        // cout<<"conf is:"<<endl;
        // for (int j = 0; j < node_num; ++ j) cout<<conf[j]<<" ";
        // cout<<endl;
        // printf("WORLD_RANK is: %d, SUB_RANK is: %d\nlocal_Hs is: %lf, Hs is: %lf\n", wd_rank, sub_rank, local_Hs, Hs);
        
    }
    MPI_Reduce(&local_Em, &Em, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ed = clock();
    if (wd_rank == 0)
    {
        Em /= conf_num;
        cout<<"The Ferromagnetism is: "<<Em<<endl;
        cout<<"time = "<<ed - st<<" s"<<endl;
    }

    MPI_Comm_free(&sub_comm);
    free(conf);
    MPI_Finalize();
    return 0;
}