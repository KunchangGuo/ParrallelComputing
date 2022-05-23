// 虚拟机下编译选项：
// aarch64-linux-gnu-g++ -pthread -fopenmp -static -o test test.cpp
// 鲲鹏服务器下编译选项：
// g++ -g -march=native -pthread -fopenmp -o test test.cpp
#include <iostream>
#include <arm_neon.h>
#include <stdlib.h>
#include <chrono>
#include <ratio>
#include <math.h>
#include <pthread.h>
#include <omp.h>
#include <unistd.h>
using namespace std;

//===线程数定义======================================================================================================================
#define NUM_THREADS 4

//===系数矩阵相关======================================================================================================================
#define MAX_N 2048
float A[MAX_N][MAX_N] = {0.0};     //在堆区申请矩阵
float A_BAC[MAX_N][MAX_N] = {0.0}; //矩阵的备份

//===函数执行方式======================================================================================================================
#define SERIAL_FUNC 0              //串行
#define DYNAMIC_FUNC 1             // pthread动态
#define STATIC_FUNC 2              // pthread减法部分静态行划分
#define STATIC_FUNC_COLUMN_MODE1 3 // pthread减法部分静态列划分(chunksize = 1)
#define STATIC_FUNC_COLUMN_MODE2 4 // pthread减法部分静态列划分(chunksize = n/num_threads)
#define SIMD_FUNC 5                // NEON
#define STATIC_FUNC_SIMD 6         // pthread减法部分静态行划分(SIMD)
#define OMP_DEFAULT_FUNC 7         // OMP默认
#define OMP_LOOP_FUNC 8            // OMP循环调度方式
#define OMP_DYNAMIC_FUNC 9         // OMP动态调度方式
#define OMP_LOOP_NEON_FUNC 10           // OMP循环调度方式+NEON
#define OMP_DYNAMIC_NEON_FUNC 11        // OMP动态调度方式+NEON

//===线程函数======================================================================================================================
void *dynamicThreadFunc(void *parm);               //动态线程函数声明
void *staticThreadFunc(void *parm);                //静态线程函数声明，消去按行划分
void *staticThreadFunc_onColumn_mode1(void *parm); //静态线程函数声明：按列划分的第一种方式：跳跃式
void *staticThreadFunc_onColumn_mode2(void *parm); //静态线程函数声明：按列划分的第二种方式：连续式
void *staticThreadFunc_SIMD(void *parm);           //静态线程函数声明：NEON版本

//===矩阵相关函数======================================================================================================================
//矩阵深拷贝
void matrixDeepCopy(int n, float b[][MAX_N], float a[][MAX_N])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            b[i][j] = a[i][j];
        }
    }
}

//矩阵初始化
void m_reset(int n, float a[][MAX_N])
{
    srand((unsigned)time(0));
    for (int i = 0; i < n; i++)
    {
        int t = pow(-1, rand() % 2 + 1);
        a[i][i] = t * rand() / float(RAND_MAX); //元素取值范围为[-1,1]之间的随机数
    }
    for (int i = 0; i < n * n; i++)
    {
        int i1 = rand() % n;
        int i2 = rand() % n;
        for (int j = 0; j < n; j++)
        {
            a[i1][j] += a[i2][j]; //随机线性组合
        }
    }
}

//矩阵显示
void printMatrix(int n, float a[][MAX_N])
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

//===串行高斯消去算法======================================================================================================================
void gaussEliminationSerial(int n, float a[][MAX_N])
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

//===线程参数定义======================================================================================================================
typedef struct
{
    int threadID; //线程ID
    int n;        //矩阵规模
    void *a;      //矩阵的空指针
} threadParam_t;

//===按行划分的动态线程消去======================================================================================================================
pthread_mutex_t mutex_dynamic;
pthread_barrier_t barrier_division_dynamic;
pthread_barrier_t barrier_elimination_dynamic;
int dynamic_row = 1;

//动态线程函数
void *dynamicThreadFunc(void *param)
{
    threadParam_t *p = (threadParam_t *)param;
    int threadID = p->threadID;                            //获取线程ID
    int n = p->n;                                          //获取矩阵规模
    float(*a)[MAX_N] = static_cast<float(*)[MAX_N]>(p->a); //获取矩阵
    int i = 1;                                             //每个线程的任务行,首先初始化为1

    for (int k = 0; k < n; k++)
    {
        if (threadID == 0) //线程0做除法
        {
            for (int j = k + 1; j < n; j++)
            {
                a[k][j] /= a[k][k];
            }
            a[k][k] = 1.0;
        }
        pthread_barrier_wait(&barrier_division_dynamic); //做完除法，所有线程同步

        while (true) //所有线程从row中取值，并递增，需要用互斥量保护
        {
            pthread_mutex_lock(&mutex_dynamic);
            i = dynamic_row;
            dynamic_row++;
            if (i >= n)
            {
                dynamic_row = k + 2;
                i = k + 2;
                pthread_mutex_unlock(&mutex_dynamic);
                sched_yield();
                break;
            }
            sched_yield();
            pthread_mutex_unlock(&mutex_dynamic);
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            a[i][k] = 0;
            i++;
        }
        pthread_barrier_wait(&barrier_elimination_dynamic); //所有线程从row中取值，并递增，需要用互斥量保护
    }
    pthread_exit(NULL);
}

//动态线程高斯消去
void gaussEliminationDynamic(int n, float a[][MAX_N])
{
    pthread_barrier_init(&barrier_division_dynamic, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_elimination_dynamic, NULL, NUM_THREADS);
    pthread_mutex_init(&mutex_dynamic, NULL);

    //创建线程并传递参数
    pthread_t handle[NUM_THREADS];          //线程句柄
    threadParam_t threadParam[NUM_THREADS]; //线程参数
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_create(&handle[threadID], NULL, dynamicThreadFunc, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    pthread_barrier_destroy(&barrier_division_dynamic);
    pthread_barrier_destroy(&barrier_elimination_dynamic);
    pthread_mutex_destroy(&mutex_dynamic);
}

//===按行划分的静态线程消去======================================================================================================================
pthread_barrier_t barrier_division_static;
pthread_barrier_t barrier_elimination_static;

//静态线程函数
void *staticThreadFunc(void *param)
{
    threadParam_t *p = (threadParam_t *)param;
    int threadID = p->threadID;                            //获取线程ID
    int n = p->n;                                          //获取矩阵规模
    float(*a)[MAX_N] = static_cast<float(*)[MAX_N]>(p->a); //获取矩阵

    for (int k = 0; k < n; k++)
    {
        if (threadID == 0)
        {
            for (int j = k + 1; j < n; j++)
            {
                a[k][j] /= a[k][k];
            }
            a[k][k] = 1.0;
        }

        pthread_barrier_wait(&barrier_division_static);

        for (int i = k + threadID + 1; i < n; i += NUM_THREADS)
        {
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            a[i][k] = 0;
        }

        pthread_barrier_wait(&barrier_elimination_static);
    }
    pthread_exit(NULL);
}

//静态线程高斯消去
void gaussEliminationStatic(int n, float a[][MAX_N])
{
    pthread_barrier_init(&barrier_division_static, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_elimination_static, NULL, NUM_THREADS);

    //创建线程并传递参数
    pthread_t handle[NUM_THREADS];          //线程句柄
    threadParam_t threadParam[NUM_THREADS]; //线程参数
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_create(&handle[threadID], NULL, staticThreadFunc, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    pthread_barrier_destroy(&barrier_division_static);
    pthread_barrier_destroy(&barrier_elimination_static);
}

//===按列划分（跳跃式）的静态线程消去======================================================================================================================
pthread_barrier_t barrier_division_static_onColumn_mode1;
pthread_barrier_t barrier_elimination_static_onColumn_mode1;

//静态线程函数,实现除法和消去均按列跳跃划分
void *staticThreadFunc_onColumn_mode1(void *param)
{
    threadParam_t *p = (threadParam_t *)param;
    int threadID = p->threadID;                            //获取线程ID
    int n = p->n;                                          //获取矩阵规模
    float(*a)[MAX_N] = static_cast<float(*)[MAX_N]>(p->a); //获取矩阵

    for (int k = 0; k < n; k++)
    {
        //计算当前线程需要做除法的任务区间
        int division_start = threadID + k + 1;
        for (int j = division_start; j < n; j += NUM_THREADS)
        {
            a[k][j] /= a[k][k];
        }

        //所有线程同步
        pthread_barrier_wait(&barrier_division_static_onColumn_mode1);

        //跳跃式分配任务
        a[k][k] = 1.0;
        for (int i = k + 1; i < n; i++)
        {
            int elimination_start = threadID + k + 1;
            for (int j = elimination_start; j < n; j += NUM_THREADS)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            //做完消去，所有线程同步
            pthread_barrier_wait(&barrier_elimination_static_onColumn_mode1);
            a[i][k] = 0;
        }
    }
    pthread_exit(NULL);
}

//静态线程高斯消去
void gaussEliminationStatic_onColumn_mode1(int n, float a[][MAX_N])
{
    pthread_barrier_init(&barrier_division_static_onColumn_mode1, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_elimination_static_onColumn_mode1, NULL, NUM_THREADS);

    //创建线程并传递参数
    pthread_t handle[NUM_THREADS];          //线程句柄
    threadParam_t threadParam[NUM_THREADS]; //线程参数
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_create(&handle[threadID], NULL, staticThreadFunc_onColumn_mode1, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    pthread_barrier_destroy(&barrier_division_static_onColumn_mode1);
    pthread_barrier_destroy(&barrier_elimination_static_onColumn_mode1);
}

//===按列划分（连续式）的静态线程消去======================================================================================================================
pthread_barrier_t barrier_division_static_onColumn_mode2;
pthread_barrier_t barrier_elimination_static_onColumn_mode2;

//静态线程函数,实现除法和消去均按列连续划分
void *staticThreadFunc_onColumn_mode2(void *param)
{
    threadParam_t *p = (threadParam_t *)param;
    int threadID = p->threadID;                            //获取线程ID
    int n = p->n;                                          //获取矩阵规模
    float(*a)[MAX_N] = static_cast<float(*)[MAX_N]>(p->a); //获取矩阵

    for (int k = 0; k < n; k++)
    {
        //计算当前线程需要做除法的任务区间
        int amount_d = n - k - 1;                                                                      //需要做除法的总任务量
        int extraTask_d = amount_d % NUM_THREADS;                                                      //额外的任务量
        int h_d = amount_d / NUM_THREADS;                                                              //步长
        int my_start_d = threadID * h_d + k + 1 + ((threadID < extraTask_d) ? threadID : extraTask_d); //当前线程的开始位置
        int my_end_d = my_start_d + h_d + ((threadID < extraTask_d) ? 1 : 0);                          //当前线程的结束位置
        for (int j = my_start_d; j < my_end_d && j < n; j++)
        {
            a[k][j] /= a[k][k];
        }

        //所有线程同步
        pthread_barrier_wait(&barrier_division_static_onColumn_mode2);

        //以线程数为步长，划分任务
        a[k][k] = 1.0;
        for (int i = k + 1; i < n; i++)
        {
            int amount_e = n - k - 1;                                                                      //需要做除法的总任务量
            int extraTask_e = amount_e % NUM_THREADS;                                                      //额外的任务量
            int h_e = amount_e / NUM_THREADS;                                                              //步长
            int my_start_e = threadID * h_e + k + 1 + ((threadID < extraTask_e) ? threadID : extraTask_e); //当前线程的开始位置
            int my_end_e = my_start_e + h_e + ((threadID < extraTask_e) ? 1 : 0);                          //当前线程的结束位置
            for (int j = my_start_e; j < my_end_e && j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }

            //做完消去，所有线程同步
            pthread_barrier_wait(&barrier_elimination_static_onColumn_mode2);
            a[i][k] = 0;
        }
    }
    pthread_exit(NULL);
}

//静态线程高斯消去
void gaussEliminationStatic_onColumn_mode2(int n, float a[][MAX_N])
{
    pthread_barrier_init(&barrier_division_static_onColumn_mode2, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_elimination_static_onColumn_mode2, NULL, NUM_THREADS);

    //创建线程
    pthread_t handle[NUM_THREADS];          //线程句柄
    threadParam_t threadParam[NUM_THREADS]; //线程参数
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_create(&handle[threadID], NULL, staticThreadFunc_onColumn_mode2, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    pthread_barrier_destroy(&barrier_division_static_onColumn_mode2);
    pthread_barrier_destroy(&barrier_elimination_static_onColumn_mode2);
}

//===NEON并行化高斯消去算法======================================================================================================================
void gaussEliminationSIMD(int n, float a[][MAX_N])
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

//===NEON、Pthread并行化高斯消去算法======================================================================================================================
pthread_barrier_t barrier_division_static_SIMD;
pthread_barrier_t barrier_elimination_static_SIMD;

// SIMD与Pthread结合的线程函数
void *staticThreadFunc_SIMD(void *param)
{
    threadParam_t *p = (threadParam_t *)param;
    int threadID = p->threadID;                            //获取线程ID
    int n = p->n;                                          //获取矩阵规模
    float(*a)[MAX_N] = static_cast<float(*)[MAX_N]>(p->a); //获取矩阵
    int offset;                                            //处理剩下元素

    float32x4_t vt;   // 保存a[k][k]
    float32x4_t va;   // 保存a[k][j]~a[k][j+3]
    float32x4_t vaik; // 保存a[i][k]
    float32x4_t vakj; // 保存a[k][j]~a[k][j+3]
    float32x4_t vaij; // 保存a[i][j]~a[i][j+3]
    float32x4_t vx;   // 保存(a[k][j]~a[k][j+3])*a[i][k]的值

    for (int k = 0; k < n; k++)
    {
        if (threadID == 0)
        {
            vt = vld1q_dup_f32(a[k] + k); // 将a[k][k]存到vt的四个通道里
            for (int j = k + 1; j + 4 <= n; j += 4)
            {
                va = vld1q_f32(a[k] + j); // 把a[k][j]~a[k][j+3]存到va的四个通道里
                va = vdivq_f32(va, vt);   // a[k][j] /= a[k][k]
                vst1q_f32(a[k] + j, va);  // 把运算结果写回到内存中
            }
            offset = (n - k - 1) % 4;
            for (int j = n - offset; j < n; j++)
            {
                a[k][j] /= a[k][k]; // 未处理的元素
            }
            a[k][k] = 1.0;
        }

        //所有线程同步
        pthread_barrier_wait(&barrier_division_static_SIMD);

        //以线程数为步长，划分任务
        for (int i = k + 1 + threadID; i < n; i += NUM_THREADS)
        {
            vaik = vld1q_f32(a[i] + k); // 将a[i][k]存到vaik的四个通道里
            for (int j = k + 1; j + 4 <= n; j += 4)
            {
                vakj = vld1q_f32(a[k] + j); // 保存a[k][j]~a[k][j+3]
                vaij = vld1q_f32(a[i] + j); // 保存a[i][j]~a[i][j+3]
                vx = vmulq_f32(vakj, vaik); // 保存(a[k][j]~a[k][j+3])*a[i][k]
                vaij = vsubq_f32(vaij, vx); // 保存(a[i][j]~a[i][j+3])-(a[k][j]~a[k][j+3])*a[i][k]
                vst1q_f32(a[i] + j, vaij);  // 将a[i][j]~a[i][j+3]写回到内存中
            }
            offset = (n - k - 1) % 4;
            for (int j = n - offset; j < n; j++)
            {
                a[i][j] -= a[k][j] * a[i][k]; // 未处理的元素
            }
            a[i][k] = 0;
        }

        //做完消去，所有线程同步
        pthread_barrier_wait(&barrier_elimination_static_SIMD);
    }
    pthread_exit(NULL);
}

// NEON与Pthread结合的高斯消去
void gaussEliminationSIMD_Pthread(int n, float a[][MAX_N])
{
    pthread_barrier_init(&barrier_division_static_SIMD, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_elimination_static_SIMD, NULL, NUM_THREADS);

    //创建线程
    pthread_t handle[NUM_THREADS];          //线程句柄
    threadParam_t threadParam[NUM_THREADS]; //线程参数
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_create(&handle[threadID], NULL, staticThreadFunc_SIMD, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < NUM_THREADS; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    pthread_barrier_destroy(&barrier_division_static_SIMD);
    pthread_barrier_destroy(&barrier_elimination_static_SIMD);
}

//===SIMD、Pthread并行化高斯消去算法======================================================================================================================

//===OpenMP高斯消去算法======================================================================================================================
void gaussEliminationOpenMP(int n, float a[][MAX_N])
{
    int i, j, k;
#pragma omp parallel num_threads(NUM_THREADS) default(none) private(i, j, k) shared(a, n)
    for (k = 0; k < n; k++)
    {
#pragma omp single
        for (j = k + 1; j < n; j++)
        {
            a[k][j] /= a[k][k];
        }
        a[k][k] = 1.0;

#pragma omp for
        for (i = k + 1; i < n; i++)
        {
            for (j = k + 1; j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            a[i][k] = 0;
        }
    }
}

//===OpenMP高斯消去，循环划分=====================================================================================================================
void gaussEliminationOpenMPLoop(int n, float a[][MAX_N])
{
    int i, j, k;
#pragma omp parallel num_threads(NUM_THREADS) default(none) private(i, j, k) shared(a, n)
    for (k = 0; k < n; k++)
    {
#pragma omp single
        for (j = k + 1; j < n; j++)
        {
            a[k][j] /= a[k][k];
        }
        a[k][k] = 1.0;

#pragma omp for schedule(static,1)
        for (i = k + 1; i < n; i++)
        {
            for (j = k + 1; j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            a[i][k] = 0;
        }
    }
}

//===OpenMP高斯消去，动态调度=====================================================================================================================
void gaussEliminationOpenMPDynamic(int n, float a[][MAX_N])
{
    int i, j, k;
#pragma omp parallel num_threads(NUM_THREADS) default(none) private(i, j, k) shared(a, n)
    for (k = 0; k < n; k++)
    {
#pragma omp single
        for (j = k + 1; j < n; j++)
        {
            a[k][j] /= a[k][k];
        }
        a[k][k] = 1.0;

#pragma omp for schedule(dynamic,1)
        for (i = k + 1; i < n; i++)
        {
            for (j = k + 1; j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            a[i][k] = 0;
        }
    }
}

//===NEON+OpenMP循环并行化高斯消去算法======================================================================================================================
void gaussElimination_OpenMP_NEON_Loop(int n, float a[][MAX_N])
{
    int k;
    int i;
    int j;
    int offset;
    
    #pragma omp parallel num_threads(NUM_THREADS) default(none) private(i, j, k, offset) shared(a, n)
    for (k = 0; k < n; k++)
    {
        #pragma omp single
        {
            float32x4_t vt = vld1q_dup_f32(a[k] + k); // 将a[k][k]存到vt的四个通道里
            for (j = k + 1; j + 4 <= n; j += 4)
            {
                float32x4_t va = vld1q_f32(a[k] + j); // 把a[k][j]~a[k][j+3]存到va的四个通道里
                va = vdivq_f32(va, vt);   // a[k][j] /= a[k][k]
                vst1q_f32(a[k] + j, va);  // 把运算结果写回到内存中
            }
            offset = (n - k - 1) % 4;
            for (j = n - offset; j < n; j++)
            {
                a[k][j] /= a[k][k]; // 未处理的元素
            }
            a[k][k] = 1.0;
        }
        
        #pragma omp for schedule(static,1)
        for (i = k + 1; i < n; i++)
        {
            float32x4_t vaik = vld1q_f32(a[i] + k); // 将a[i][k]存到vaik的四个通道里
            for (j = k + 1; j + 4 <= n; j += 4)
            {
                float32x4_t vakj = vld1q_f32(a[k] + j); // 保存a[k][j]~a[k][j+3]
                float32x4_t vaij = vld1q_f32(a[i] + j); // 保存a[i][j]~a[i][j+3]
                float32x4_t vx = vmulq_f32(vakj, vaik); // 保存(a[k][j]~a[k][j+3])*a[i][k]
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

//===NEON+OpenMP循环并行化高斯消去算法======================================================================================================================
void gaussElimination_OpenMP_NEON_Dynamic(int n, float a[][MAX_N])
{
    int k;
    int i;
    int j;
    int offset;
    
    #pragma omp parallel num_threads(NUM_THREADS) default(none) private(i, j, k, offset) shared(a, n)
    for (k = 0; k < n; k++)
    {
        #pragma omp single
        {
            float32x4_t vt = vld1q_dup_f32(a[k] + k); // 将a[k][k]存到vt的四个通道里
            for (j = k + 1; j + 4 <= n; j += 4)
            {
                float32x4_t va = vld1q_f32(a[k] + j); // 把a[k][j]~a[k][j+3]存到va的四个通道里
                va = vdivq_f32(va, vt);   // a[k][j] /= a[k][k]
                vst1q_f32(a[k] + j, va);  // 把运算结果写回到内存中
            }
            offset = (n - k - 1) % 4;
            for (j = n - offset; j < n; j++)
            {
                a[k][j] /= a[k][k]; // 未处理的元素
            }
            a[k][k] = 1.0;
        }
        
        #pragma omp for schedule(dynamic,1)
        for (i = k + 1; i < n; i++)
        {
            float32x4_t vaik = vld1q_f32(a[i] + k); // 将a[i][k]存到vaik的四个通道里
            for (j = k + 1; j + 4 <= n; j += 4)
            {
                float32x4_t vakj = vld1q_f32(a[k] + j); // 保存a[k][j]~a[k][j+3]
                float32x4_t vaij = vld1q_f32(a[i] + j); // 保存a[i][j]~a[i][j+3]
                float32x4_t vx = vmulq_f32(vakj, vaik); // 保存(a[k][j]~a[k][j+3])*a[i][k]
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

//===计时函数======================================================================================================================
double getTime(int n, float a[][MAX_N], int mode)
{
    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();
    switch (mode)
    {
    case SERIAL_FUNC:
        gaussEliminationSerial(n, a);
        break;
    case DYNAMIC_FUNC:
        gaussEliminationDynamic(n, a);
        break;
    case STATIC_FUNC:
        gaussEliminationStatic(n, a);
        break;
    case STATIC_FUNC_COLUMN_MODE1:
        gaussEliminationStatic_onColumn_mode1(n, a);
        break;
    case STATIC_FUNC_COLUMN_MODE2:
        gaussEliminationStatic_onColumn_mode2(n, a);
        break;
    case SIMD_FUNC:
        gaussEliminationSIMD(n, a);
        break;
    case STATIC_FUNC_SIMD:
        gaussEliminationSIMD_Pthread(n, a);
        break;
    case OMP_DEFAULT_FUNC:
        gaussEliminationOpenMP(n, a);
        break;
    case OMP_LOOP_FUNC:
        gaussEliminationOpenMPLoop(n, a);
        break;
    case OMP_DYNAMIC_FUNC:
        gaussEliminationOpenMPDynamic(n, a);
        break;
    case OMP_LOOP_NEON_FUNC:
        gaussElimination_OpenMP_NEON_Loop(n, a);
        break;
    case OMP_DYNAMIC_NEON_FUNC:
        gaussElimination_OpenMP_NEON_Dynamic(n, a);
        break;
    default:
        break;
    }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    return time_span.count();
}

//===主函数======================================================================================================================
int main()
{
    int n = 1024;
    m_reset(n, A);
    matrixDeepCopy(n, A_BAC, A);

    cout << "Serial:         " << getTime(n, A, SERIAL_FUNC) << endl;
    // printMatrix(n, A);

    matrixDeepCopy(n, A, A_BAC);
    cout << "OpenMP:         " << getTime(n, A, OMP_DEFAULT_FUNC) << endl;
    // printMatrix(n, A);

    matrixDeepCopy(n, A, A_BAC);
    cout << "Pthread:        " << getTime(n, A, STATIC_FUNC) << endl;
    // printMatrix(n, A);

    matrixDeepCopy(n, A, A_BAC);
    cout << "NEON:           " << getTime(n, A, SIMD_FUNC) << endl;
    // printMatrix(n, A);

    matrixDeepCopy(n, A, A_BAC);
    cout << "NEON & Pthread: " << getTime(n, A, STATIC_FUNC_SIMD) << endl;
    // printMatrix(n, A);

    matrixDeepCopy(n, A, A_BAC);
    cout << "NEON & OpenMPs: " << getTime(n, A, OMP_LOOP_NEON_FUNC) << endl;

    matrixDeepCopy(n, A, A_BAC);
    cout << "NEON & OpenMPd: " << getTime(n, A, OMP_DYNAMIC_NEON_FUNC) << endl;

    return 0;
}