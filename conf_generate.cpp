//多进程并行写多个文件
// #include<mpi.h>
// #include<stdio.h>
// #include<stdlib.h>
// #include<cstring>
// #include<iostream>
// #include<cmath>
// using namespace std;

// void get_binary(char *vt, int x, int node_num)
// {
//     int i = node_num - 1;
//     while (x)
//     {
//         vt[i --] = (x & 1) + '0';
//         x >>= 1;
//     }
//     while (i >= 0) vt[i --] = '0';
// }

// int main()
// {
//     int comm_rk, comm_sz;
//     MPI_Init(NULL, NULL);
//     MPI_Comm_rank(MPI_COMM_WORLD, &comm_rk);
//     MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

//     int node_num = 2;
//     int conf_num = pow(2, node_num);
//     char *buf = (char *)calloc(node_num, sizeof(char));

//     int local_n = conf_num / comm_sz;
//     int local_st = comm_rk * local_n;
//     int local_ed = local_st + local_n - 1;

//     MPI_Status status;
//     MPI_File fh;

//     for (int i = local_st; i <= local_ed; ++ i)
//     {
//         get_binary(buf, i, node_num);
//         string g = "/home/ntyz/work_projects/MPI_exercises/confs_data/" + to_string(i) + ".txt";
//         MPI_File_open(MPI_COMM_SELF, g.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        
//         MPI_File_write(fh, buf, strlen(buf), MPI_CHAR, &status);
//     }

//     MPI_File_close(&fh);
//     free(buf);
//     MPI_Finalize();
//     return 0;  
// } 


// 多个进程写单个文件方法
#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<cstring>
#include<iostream>
#include<cmath>
using namespace std;

void get_binary(char *vt, int x, int node_num)
{
    int i = node_num - 1;
    while (x)
    {
        vt[i --] = (x & 1) + '0';
        x >>= 1;
    }
    while (i >= 0) vt[i --] = '0';
}

int main()
{
    int comm_rk, comm_sz;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rk);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    int node_num = 2;
    int conf_num = pow(2, node_num);

    int local_n = conf_num / comm_sz;
    int local_st = comm_rk * local_n;
    int local_ed = local_st + local_n - 1;

    char *buf = (char *)calloc(node_num * local_n, sizeof(char));

    // MPI_Offset offset = comm_rk * 256;
    MPI_Request request;
    MPI_Status status;
    MPI_File fh;
    MPI_Offset offset = comm_rk * local_n;
    string g = "/home/ntyz/work_projects/MPI_exercises/confs_data/confs_file.txt";
    MPI_File_open(MPI_COMM_SELF, g.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    
    for (int i = local_st; i <= local_ed; ++ i)
    {
        get_binary(buf, i, node_num);
        MPI_File_iwrite_at(fh, (offset + (i%local_n)) * node_num * sizeof(char), buf, node_num, MPI_CHAR, &request);
	    MPI_Wait(&request, &status);
    }

    // cout<<"the length of buf is "<<strlen(buf)<<endl;
    // cout<<"the base offset is "<<offset<<endl;
    // char *rbuf = (char *)calloc(node_num, sizeof(char));
    
    // MPI_Offset poffset = offset * strlen(buf);
    MPI_File_iread_at(fh, offset * node_num * sizeof(char), buf, node_num * local_n, MPI_CHAR, &request);  // strlen(buf) * local_n
	MPI_Wait(&request, &status);

    // cout<<"the process offset is "<<poffset<<endl;
    // cout<<"the read length of rbuf is "<<strlen(buf) * local_n<<endl;
    // cout<<"the length of rbuf is "<<strlen(rbuf)<<endl;
    printf("rankID = %d, content = \'%s\'\n", comm_rk, buf);

    MPI_File_close(&fh);
    free(buf); 
    MPI_Finalize();
    return 0;  
} 