#include <cstdio>
#include<iostream>
#include<string>
#include<cstring>
#include<fstream>
#include<cmath>
#include<cstdlib>
using namespace std;

// test send value to c++ file by arguements 
int main(int argc, char **argv)
{
    int sum = 0;
    for(int i = 1; i < argc; ++ i)
    {
        int num_i = atoi(argv[i]);
        sum += num_i;
    } 
    cout<<"result is: "<<sum<<endl;
    return 0;
}
// 第一个参数是当前运行的文件名，所以索引从 1 开始。 运行方式为 ./a.exe 3 2
// 编译运行方法
// g++ -o C++_args C++_args.cpp 
// ./C++_args 1 2 3 4 5