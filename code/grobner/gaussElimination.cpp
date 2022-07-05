#include <iostream>
#include <arm_neon.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
using namespace std;

const int maxN = 640; // 系数矩阵最大规模

// 矩阵深拷贝:a拷贝给b
void matrixDeepCopy(int n, float b[][maxN], float a[][maxN])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            b[i][j] = a[i][j];
        }
    }
}

// 串行高斯消去
void gaussEliminationOrdinary(int n, float a[][maxN])
{
    for (int k = 0; k < n; k++)
    {
        for (int j = k + 1; j < n; j++)
        {
            a[k][j] /= a[k][k];
        }
        a[k][k] = 1.0;
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            a[i][k] = 0;
        }
    }
}

// SIMD化的高斯消去
void gaussEliminationOptimisedUnmatched(int n, float a[][maxN])
{
    float32x4_t vt;   // 保存a[k][k]
    float32x4_t va;   // 保存a[k][j]~a[k][j+3]
    float32x4_t vaik; // 保存a[i][k]
    float32x4_t vakj; // 保存a[k][j]~a[k][j+3]
    float32x4_t vaij; // 保存a[i][j]~a[i][j+3]
    float32x4_t vx;   // 保存(a[k][j]~a[k][j+3])*a[i][k]的值

    int k;
    int i;
    int j;
    int offset;
    for (k = 0; k < n; k++)
    {
        vt = vld1q_dup_f32(a[k] + k); // 将a[k][k]存到vt的四个通道里
        for (j = k + 1; j + 4 <= n; j += 4)
        {
            va = vld1q_f32(a[k] + j); // 把a[k][j]~a[k][j+3]存到va的四个通道里
            va = vdivq_f32(va, vt);   // a[k][j] /= a[k][k]
            vst1q_f32(a[k] + j, va);  // 把运算结果写回到内存中
        }
        offset = (n - k - 1) % 4;
        for (j = n - offset; j < n; j++)
        {
            a[k][j] /= a[k][k]; // 未处理的元素
        }
        a[k][k] = 1.0;
        for (i = k + 1; i < n; i++)
        {
            vaik = vld1q_f32(a[i] + k); // 将a[i][k]存到vaik的四个通道里
            for (j = k + 1; j + 4 <= n; j += 4)
            {
                vakj = vld1q_f32(a[k] + j); // 保存a[k][j]~a[k][j+3]
                vaij = vld1q_f32(a[i] + j); // 保存a[i][j]~a[i][j+3]
                vx = vmulq_f32(vakj, vaik); // 保存(a[k][j]~a[k][j+3])*a[i][k]
                vaij = vsubq_f32(vaij, vx); // 保存(a[i][j]~a[i][j+3])-(a[k][j]~a[k][j+3])*a[i][k]
                vst1q_f32(a[i] + j, vaij);  // 将a[i][j]~a[i][j+3]写回到内存中
            }
            for (j = n - offset; j < n; j++)
            {
                a[i][j] -= a[k][j] * a[i][k]; // 未处理的元素
            }
            a[i][k] = 0;
        }
    }
}

// 对齐的SIMD化的高斯消去
void gaussEliminationOptimisedMatched(int n, float a[][maxN])
{
    float32x4_t vt;   // 保存a[k][k]
    float32x4_t va;   // 保存a[k][j]~a[k][j+3]
    float32x4_t vaik; // 保存a[i][k]
    float32x4_t vakj; // 保存a[k][j]~a[k][j+3]
    float32x4_t vaij; // 保存a[i][j]~a[i][j+3]
    float32x4_t vx;   // 保存(a[k][j]~a[k][j+3])*a[i][k]的值

    int k;
    int i;
    int j;
    int start;
    int leftover = n % 4;
    int end = n - leftover;
    int handle;
    for (k = 0; k < n; k++)
    {
        vt = vld1q_dup_f32(a[k] + k); // 将a[k][k]存到vt的四个通道里
        start = k + 1;
        handle = start % 4;
        // 掐头
        switch (handle)
        {
        case 0:
            break; // 开始位置已经对齐
        case 1:
            a[k][start] /= a[k][k];
            a[k][start + 1] /= a[k][k];
            a[k][start + 2] /= a[k][k];
            start += 3;
            break;
        case 2:
            a[k][start] /= a[k][k];
            a[k][start + 1] /= a[k][k];
            start += 2;
            break;
        case 3:
            a[k][start] /= a[k][k];
            start += 1;
            break;
        }
        for (j = start; j + 4 <= end; j += 4)
        {
            va = vld1q_f32(a[k] + j); // 把a[k][j]~a[k][j+3]存到va的四个通道里
            va = vdivq_f32(va, vt);   // a[k][j] /= a[k][k]
            vst1q_f32(a[k] + j, va);  // 把运算结果写回到内存中
        }
        // 去尾巴
        for (j = end; j < n; j++)
        {
            a[k][j] /= a[k][k]; // 尾巴部分
        }
        a[k][k] = 1.0;
        for (i = k + 1; i < n; i++)
        {
            vaik = vld1q_f32(a[i] + k); // 将a[i][k]存到vaik的四个通道里
            start = k + 1;
            handle = start % 4;
            // 掐头
            switch (handle)
            {
            case 0:
                break; // 开始位置已经对齐
            case 1:
                a[i][start] -= a[i][k] * a[k][start];
                a[i][start + 1] -= a[i][k] * a[k][start + 1];
                a[i][start + 2] -= a[i][k] * a[k][start + 2];
                start += 3;
                break;
            case 2:
                a[i][start] -= a[i][k] * a[k][start];
                a[i][start + 1] -= a[i][k] * a[k][start + 1];
                start += 2;
                break;
            case 3:
                a[i][start] -= a[i][k] * a[k][start];
                start += 1;
                break;
            }
            for (j = start; j + 4 <= end; j += 4)
            {
                vakj = vld1q_f32(a[k] + j); // 保存a[k][j]~a[k][j+3]
                vaij = vld1q_f32(a[i] + j); // 保存a[i][j]~a[i][j+3]
                vx = vmulq_f32(vakj, vaik); // 保存(a[k][j]~a[k][j+3])*a[i][k]
                vaij = vsubq_f32(vaij, vx); // 保存(a[i][j]~a[i][j+3])-(a[k][j]~a[k][j+3])*a[i][k]
                vst1q_f32(a[i] + j, vaij);  // 将a[i][j]~a[i][j+3]写回到内存中
            }
            //  去尾
            for (j = end; j < n; j++)
            {
                a[i][j] -= a[k][j] * a[i][k]; // 未处理的元素
            }
            a[i][k] = 0;
        }
    }
}

void printMatrix(int n, float a[][maxN])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << a[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}

void m_reset(int n, float a[][maxN])
{
    srand((unsigned)time(0));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            a[i][j] = 0;
        }
        a[i][i] = 1.0;
        for (int j = i + 1; j < n; j++)
        {
            int t = pow(-1, rand() % 2 + 1);
            a[i][j] = t * rand() / float(RAND_MAX);
        }
    }
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                a[i][j] += a[k][j];
            }
        }
    }
}

int main()
{
    int n;           // 矩阵为n*n
    int times = 3;   // 重复测量次数为3
    int amount = 16; // 测量的矩阵规模有16组

    /*int sizeN[] = {2,3,4,5,
        6,7,8,12,
        16,20,24,28,
        32,46,47,48,
        49,64,96,128,
        196,256,512,1024};*/
    int sizeN[] = {5, 72, 80, 88, 96, 104, 112, 120, 125, 126, 127, 128, 129, 130, 131, 132};
    for (int i = 0; i < amount; i++)
    {
        n = sizeN[i];

        float a[n][maxN] = {0};
        float b[n][maxN] = {0};
        float c[n][maxN] = {0};
        m_reset(n, a);           // 生成一个a
        matrixDeepCopy(n, b, a); // b为a的复制
        matrixDeepCopy(n, c, a); // b为a的复制
        // printMatrix(n,a);

        // Linux下高精度计时,对每个规模都重复5次来计量以减少计算误差
        time_t dsec_a = 0, dsec_b = 0, dsec_c = 0;
        long long dnsec_a = 0, dnsec_b = 0, dnsec_c = 0;

        for (int j = 0; j < times; j++)
        {
            // 测量未对齐的SIMD化的高斯消去时间
            struct timespec sts, ets;
            timespec_get(&sts, TIME_UTC);             // 测量的开始时间
            gaussEliminationOptimisedUnmatched(n, a); // 测量
            timespec_get(&ets, TIME_UTC);             // 测量的结束时间
            time_t dsec = ets.tv_sec - sts.tv_sec;
            long dnsec = ets.tv_nsec - sts.tv_nsec;
            if (dnsec < 0)
            {
                dsec--;
                dnsec += 1000000000ll;
            }
            // printMatrix(n,a);
            // printf ("%lld.%09llds\n",dsec,dnsec);
            dsec_a += dsec; // 将测量的时间进行累加
            dnsec_a += dnsec;

            // 测量未并行的高斯消去时间
            timespec_get(&sts, TIME_UTC);   // 测量的开始时间
            gaussEliminationOrdinary(n, b); // 测量
            timespec_get(&ets, TIME_UTC);   // 测量的结束时间
            dsec = ets.tv_sec - sts.tv_sec;
            dnsec = ets.tv_nsec - sts.tv_nsec;
            if (dnsec < 0)
            {
                dsec--;
                dnsec += 1000000000ll;
            }
            // printMatrix(n,a);
            // printf ("%lld.%09llds\n",dsec,dnsec);
            dsec_b += dsec; // 将测量的时间进行累加
            dnsec_b += dnsec;

            // 测量对齐的SIMD高斯消去时间
            timespec_get(&sts, TIME_UTC);           // 测量的开始时间
            gaussEliminationOptimisedMatched(n, c); // 测量
            timespec_get(&ets, TIME_UTC);           // 测量的结束时间
            dsec = ets.tv_sec - sts.tv_sec;
            dnsec = ets.tv_nsec - sts.tv_nsec;
            if (dnsec < 0)
            {
                dsec--;
                dnsec += 1000000000ll;
            }
            // printMatrix(n,a);
            // printf ("%lld.%09llds\n",dsec,dnsec);
            dsec_c += dsec; // 将测量的时间进行累加
            dnsec_c += dnsec;

            // 重新初始化a矩阵和b矩阵
            m_reset(n, a);
            matrixDeepCopy(n, b, a);
            matrixDeepCopy(n, c, a);
        }

        // 获得平均用时
        dsec_a /= times;
        dsec_b /= times;
        dnsec_a /= times;
        dnsec_b /= times;
        dsec_c /= times;
        dnsec_c /= times;

        printf("ordinary   n= %d: ordinary  = %11d.%09llds\n", n, dsec_b, dnsec_b);
        printf("optimised  n= %d: optimised = %11d.%09llds\n", n, dsec_a, dnsec_a);
        printf("matched    n= %d: optimised = %11d.%09llds\n", n, dsec_c, dnsec_c);
    }
    return 0;
}
