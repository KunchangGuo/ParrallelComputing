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
//void init()
//{
//    for(int i = 0; i < N; i++)
//    {
//        A[i] = i;
//    }
//}

//sum
void sumarize()
{
    for(int i = 0; i < N; i++)
    {
        sum += A[i];
    }
}

int main()
{
    long long head, tail, freq;

    //times
    int counter = 0;

    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    //frequency
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    do
    {
        counter++;
//        init();
        sumarize();

        //end time
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);

    }while((tail - head) * 1.0 / freq < 3);

    //display
    cout<<"Loops: "<<counter<<endl;
    cout<<"Time: "<<(tail - head) * 1000.0 / freq<<endl;
    cout<<"TimePerLoop: "<<(tail - head) * 1000.0 / (freq * counter)<<endl;

    return 0;
}
