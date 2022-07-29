#include <cstdio>
#include<iostream>
#include<string>
#include<cstring>
#include<fstream>
#include<cmath>
#include<cstdlib>
using namespace std;

// test read/write file by c++
int main() 
{ 
    ifstream myfile("/home/ntyz/work_projects/MPI_exercises/in.txt"); 
    ofstream outfile("/home/ntyz/work_projects/MPI_exercises/out.txt", ios::app); 
    string temp; 
    if (!myfile.is_open()) 
    { 
        cout << "未成功打开文件" << endl; 
    } 
    while(getline(myfile,temp)) 
    { 
        outfile << temp; 
        outfile << endl;
    } 
    myfile.close(); 
    outfile.close();
    return 0; 
} 