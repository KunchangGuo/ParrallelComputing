#include <iostream>
#include <windows.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

//scale
const int N = pow(2, 27);

double A[N];
double sum;

//initialize
void init(int n)
{
    sum = 0.0;
    for(int i = 0; i < n; i++)
    {
        A[i] = i;
    }
}

//sum
void sumarize(int n)
{
    for(int i = 0; i < n; i++)
    {
        sum += A[i];
    }
}

// 多链路式
void multiWays(int n)
{

    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < n; i += 2)
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
        return;
    }
    else
    {
        for(int i = 0; i < n / 2; i++)
        {
            A[i] += A[n-1-i];
        }
        n /= 2;
        recursion(n);
    }
}

//二重循环
void doubleLoop(int n)
{
    for(int m = n; m > 1; m /= 2)
    {
        for(int i = 0; i < m/2; i++)
        {
            A[i] = A[2*i] + A[2*i+1];
        }
    }
}


int main()
{
    long long head, tail, freq;

    //times
    int counter = 0;

    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    init(N);

    //frequency
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    do
    {
        counter++;

        //normal
//        sumarize(N);

//        multiWays(N);

//        recursion(N);

        doubleLoop(N);

        //end time
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);

    }while((tail - head) * 1.0 / freq < 3);

    //display
    cout<<"Loops: "<<counter<<endl;
    cout<<"Time: "<<(tail - head) * 1000.0 / freq<<endl;
    cout<<"TimePerLoop: "<<(tail - head) * 1000.0 / (freq * counter)<<endl;

    return 0;
}
