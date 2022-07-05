#include <emmintrin.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
using namespace std;

int main()
{
    int a[8] = {0,0,0,0,1,0,1,1};
    int b[8] = {0,0,0,0,1,1,0,1};
    __m128i v1 = _mm_load_si128((__m128i*)(a+4));
    __m128i v2 = _mm_load_si128((__m128i*)(b+4));
    // result = {0,0,0,0,0,0,1,1}
    v1 = _mm_xor_si128(v1,v2);
    _mm_store_si128((__m128i*)(a+4),v1);
    for(int i=0;i<8;i++)
    {
        cout<<a[i]<<' ';
    }
    cout<<endl;
    return 0;
}