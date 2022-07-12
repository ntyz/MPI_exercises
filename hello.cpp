#include<iostream>
#include<mpi.h>
#include<unistd.h>
using namespace std;

int main()
{
    int i = 0;
    while (i == 0)
    {
        sleep(5);
    }
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char data = 'A';
    if (rank == 0)
    {
        data = 'B';
        MPI_Send(&data, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        cout<<"rank is: "<<rank<<" data is: "<<data<<endl;
    }
    else if (rank == 1)
    {
        MPI_Recv(&data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout<<"rank is: "<<rank<<" data is: "<<data<<endl;
    }
    return 0;
}