#include <iostream>
#include <windows.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

//scale
const int N = pow(2, 21);

double A[N];
double tmp[N];
double sum;

//initialize
void init()
{
    sum = 0.0;
    for(int i = 0; i < N; i++)
    {
        A[i] = i;
        tmp[i] = 0.0;
    }
}

//sum
void summarize()
{
    for(int i = 0; i < N; i++)
    {
        sum += A[i];
    }
}

// 多链路式
void multiWays()
{

    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < N; i += 4)
    {
        sum1 += A[i];
        sum2 += A[i + 1];
    }
    sum = sum1 + sum2;
}

//递归
void recursion(int n)
{
    if(n == 1)
    {
        sum = tmp[0];
        return;
    }
    else
    {
        for(int i = 0; i < n / 2; i++)
        {
//            A[i] += A[n-1-i];
            tmp[i] = A[i] + A[n-1-i];
        }
        n /= 2;
        recursion(n);
    }
}

//二重循环
void doubleLoop(int N)
{
    for(int m = N; m > 1; m /= 2)
    {
        for(int i = 0; i < m/2; i++)
        {
//            A[i] = A[2*i] + A[2*i+1];
            tmp[i] = A[2*i] + A[2*i+1];
        }
    }
    sum = tmp[0];
}


int main()
{
    long long commonHead, commonTail;
    long long optimisedHead, optimisedTail;

    //frequency
    long long freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    //times
    int commonCounter = 0;
    int optimisedCounter = 0;

    init();

    QueryPerformanceCounter((LARGE_INTEGER*)&commonHead);
    do
    {
        commonCounter++;
        sum = 0;

        summarize();

        //end time
        QueryPerformanceCounter((LARGE_INTEGER*)&commonTail);

    }while(double(commonTail - commonHead) / freq < 10.0);


    QueryPerformanceCounter((LARGE_INTEGER*)&optimisedHead);
    do
    {
        optimisedCounter++;
        sum = 0;

        multiWays();

        //end time
        QueryPerformanceCounter((LARGE_INTEGER*)&optimisedTail);

    }while(double(optimisedTail - optimisedHead) / freq < 10.0);

    //display
    cout<<"CommonCounter: "<<commonCounter<<endl
        <<"CommonDuration: "<<(commonTail - commonHead)*1000000.0 / freq<<endl
        <<"CommonPeriod: "<<(commonTail - commonHead) * 1000000.0 / (freq * commonCounter)<<endl<<endl;
    cout<<"OptimisedCounter: "<<optimisedCounter<<endl
        <<"OptimisedDuration: "<<(optimisedTail - optimisedHead)*1000000.0 / freq<<endl
        <<"OptimisedPeriod: "<<(optimisedTail - optimisedHead) * 1000000.0 / (freq * optimisedCounter)<<endl<<endl;

    return 0;
}
