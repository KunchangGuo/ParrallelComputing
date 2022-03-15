#include <iostream>
#include <windows.h> //高精度计时
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
using namespace std;

//scale
const int N = pow(2,11);

double A[N][N], B[N], C[N];

void init()
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
             A[i][j] = i + j;
        }
        B[i] = i;
    }
}

void commonAlgorithm()
{
    //visiting by row
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            C[i] += A[j][i] * B[j];
        }
    }
}

void optimisedAlgorithm()
{
    //visiting by column
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            C[i] += A[j][i] * B[j];
        }
    }
}

int main()
{

    //start time and end time
    long long commonHead, commonTail;
    long long optimisedHead, optimisedTail;

    frequency of CPU clock
    long long freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    times
    int commonCounter = 0;
    int optimisedCounter = 0;

    measuring of common algorithm
    QueryPerformanceCounter((LARGE_INTEGER*)&commonHead);

    do
    {
        commonCounter++;
        for(int i = 0; i < N; i++)
        {
            C[i] = 0.0;
        }
        commonAlgorithm();

        QueryPerformanceCounter((LARGE_INTEGER*)&commonTail);

    }while((double)(commonTail-commonHead) / freq < 1.0);

    //measuring of optimised algorithm
    QueryPerformanceCounter((LARGE_INTEGER*)&optimisedHead);
    do
    {
        optimisedCounter++;
        for(int i = 0; i < N; i++)
        {
            C[i] = 0.0;
        }
        optimisedAlgorithm();

        QueryPerformanceCounter((LARGE_INTEGER*)&optimisedTail);

    }while(double(optimisedTail - optimisedHead) / freq < 10.0);

    //show results
    cout<<"CommonTimes: "<<commonCounter<<endl
        <<"CommonDuration: "<<(commonTail - commonHead) * 1000000.0 / freq<<" us"<<endl
        <<"CommonPeriod: "<<(commonTail - commonHead) * 1000000.0 / (freq * commonCounter)<<" us per loop"<<endl<<endl;

    cout<<"OptimisedTimes: "<<optimisedCounter<<endl
    <<"OptimisedDuration: "<<(optimisedTail - optimisedHead) * 1000000.0 / freq<<" us"<<endl
    <<"OptimisedPeriod: "<<(optimisedTail - optimisedHead) * 1000000.0 / (freq * optimisedCounter)<<" us per loop"<<endl;

    return 0;
}
