//多进程并行读写一个文件
// 多个进程写单个文件方法
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>
 
/**
 * all process read/write at the same offset
 * int MPI_File_iread_at(MPI_File fh, MPI_Offset offset, void * buf, int count, MPI_Datatype datatype, MPI_Request * request)
 *		after call, the process immediately go back to executing on any memory except the buf area
 * MPI_Wait(&request, &status) is used to check whether the data in buf is ready
 * int MPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, void * buf, int count, MPI_Datatype datatype, MPI_Request * request)
 * int MPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void * buf, int count, MPI_Datatype datatype)
 */
 
int main(int argc, char *argv[])
{
	int totalTaskNum, rankID;
 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rankID);
	MPI_Comm_size(MPI_COMM_WORLD, &totalTaskNum);
 
	char *filename = "file.txt";
	MPI_File fh;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
 
	MPI_Offset offset = rankID * 256;
	MPI_Request request;
	char buf[100];
	sprintf(buf, "rankID = %d, totalTaskNum = %d, Hello World!\n", rankID, totalTaskNum);
	MPI_File_iwrite_at(fh, offset, buf, strlen(buf), MPI_CHAR, &request);
	MPI_Wait(&request, &status);
 
 
	MPI_File_iread_at(fh, offset, buf, sizeof(buf), MPI_CHAR, &request);
	MPI_Wait(&request, &status);
	printf("rankID = %d, content = \'%s\'\n", rankID, buf);
 
	MPI_File_close(&fh);   //after open, fh has the communicator info
 
	MPI_Finalize();
	return 0;
}




//多进程并行写多个文件

// #include<mpi.h>
// #include<stdio.h>
// #include<stdlib.h>
// #include<cstring>
// #include<iostream>
// #include<cmath>
// using namespace std;

// int main()
// {
//     int comm_rk, comm_sz;
//     MPI_Init(NULL, NULL);
//     MPI_Comm_rank(MPI_COMM_WORLD, &comm_rk);
//     MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

//     MPI_Status status;
//     MPI_File fh;
// 	string test;
// 	const char* buf[100];
// 	string g = "/home/ntyz/work_projects/MPI_exercises/confs_data/" + to_string(comm_rk) + ".txt";
// 	MPI_File_open(MPI_COMM_SELF, g.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

// 	for (int i = 0; i < 100; ++ i)
// 	{
// 		test = to_string(i) + "\n";
// 		buf[i] = test.c_str();
// 		MPI_File_write(fh, buf[i], strlen(buf[i]), MPI_CHAR, &status);
// 	}
	
//     MPI_File_close(&fh);
//     MPI_Finalize();
//     return 0;  
// } 