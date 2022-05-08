#include <iostream>
#include <arm_neon.h>
#include <stdlib.h>
#include <chrono>
#include <ratio>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
using namespace std;

//===线程数定义和稀疏矩阵最大规模定义以及函数指针声明======================================================================================================================

#define THREAD_NUM 4 //线程数
#define MAX_N 520    //系数矩阵最大规模

#define SERIAL_FUNC 0
#define DYNAMIC_FUNC 1
#define STATIC_FUNC 2
#define STATIC_FUNC_COLUMN_MODE1 3
#define STATIC_FUNC_COLUMN_MODE2 4
#define SIMD_FUNC 5
#define STATIC_FUNC_SIMD 6

void *dynamicThreadFunc(void *parm);               //动态线程函数声明
void *staticThreadFunc(void *parm);                //静态线程函数声明，消去按行划分
void *staticThreadFunc_onColumn_mode1(void *parm); //静态线程函数声明：按列划分的第一种方式：跳跃式
void *staticThreadFunc_onColumn_mode2(void *parm); //静态线程函数声明：按列划分的第二种方式：连续式
void *staticThreadFunc_SIMD(void *parm);           //静态线程函数声明：SIMD版本

//===线程数定义和稀疏矩阵最大规模定义以及函数指针声明======================================================================================================================

//===有关矩阵的函数操作======================================================================================================================

//矩阵深拷贝:a拷贝给b，以相同矩阵进行运算以测量不同策略下的运算时间
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

//初始化矩阵
void m_reset(int n, float a[][MAX_N])
{
    srand((unsigned)time(0));
    //初始化上三角矩阵,元素取值范围为[-1,1]之间的随机数
    for (int i = 0; i < n; i++)
    {
        int t = pow(-1, rand() % 2 + 1);
        a[i][i] = t * rand() / float(RAND_MAX);
    }
    //随机选择两行进行线性相加，重复n^2次
    for (int i = 0; i < n * n; i++)
    {
        int i1 = rand() % n;
        int i2 = rand() % n;
        for (int j = 0; j < n; j++)
        {
            a[i1][j] += a[i2][j];
        }
    }
}

//打印矩阵
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

//===有关矩阵的函数操作======================================================================================================================

//===串行高斯消去算法======================================================================================================================

//串行高斯消去
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

//===串行高斯消去算法======================================================================================================================

//===线程数据结构定义======================================================================================================================

//线程数据结构定义
typedef struct
{
    int threadID; //线程ID
    int n;        //矩阵规模
    void *a;      //矩阵
} threadParam_t;

//===线程数据结构定义======================================================================================================================

//===按行划分的动态线程消去======================================================================================================================

//声明互斥量
pthread_mutex_t mutex_dynamic;
//声明barrier
pthread_barrier_t barrier_division_dynamic;
pthread_barrier_t barrier_elimination_dynamic;
//定义动态线程消去的行
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
        //线程0做除法，其他线程等待
        if (threadID == 0)
        {
            for (int j = k + 1; j < n; j++)
            {
                a[k][j] /= a[k][k];
            }
            a[k][k] = 1.0;
        }
        // cout << "From Thread :" << threadID << endl;
        //做完除法，所有线程同步
        pthread_barrier_wait(&barrier_division_dynamic);
        // cout << "done first step!" << endl;

        //所有线程从row中取值，并递增，需要用互斥量保护
        while (true)
        {
            //有保护地获取当前任务行
            pthread_mutex_lock(&mutex_dynamic);
            i = dynamic_row;
            // cout << "From Thread :" << threadID << " i = " << i << endl;
            dynamic_row++;
            if (i >= n)
            {
                //下一次k取k+1，要从k+2开始消去
                dynamic_row = k + 2;
                i = k + 2;
                pthread_mutex_unlock(&mutex_dynamic);
                sched_yield();
                break;
            }
            sched_yield();
            pthread_mutex_unlock(&mutex_dynamic);
            //执行消去
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            a[i][k] = 0;
            i++;
        }
        //做完消去，所有线程同步
        pthread_barrier_wait(&barrier_elimination_dynamic);
    }
    pthread_exit(NULL);
}

//动态线程高斯消去
void gaussEliminationDynamic(int n, float a[][MAX_N])
{
    //初始化barrier
    pthread_barrier_init(&barrier_division_dynamic, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_elimination_dynamic, NULL, THREAD_NUM);

    //初始化互斥量
    pthread_mutex_init(&mutex_dynamic, NULL);

    //创建线程
    pthread_t handle[THREAD_NUM];          //线程句柄
    threadParam_t threadParam[THREAD_NUM]; //线程参数
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_create(&handle[threadID], NULL, dynamicThreadFunc, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    //销毁barrier
    pthread_barrier_destroy(&barrier_division_dynamic);
    pthread_barrier_destroy(&barrier_elimination_dynamic);
    //销毁互斥量
    pthread_mutex_destroy(&mutex_dynamic);
}

//===按行划分的动态线程消去======================================================================================================================

//===按行划分的静态线程消去======================================================================================================================

//定义barrier
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
        //线程0做除法，其他线程等待
        if (threadID == 0)
        {
            for (int j = k + 1; j < n; j++)
            {
                a[k][j] /= a[k][k];
            }
            a[k][k] = 1.0;
        }

        //所有线程同步
        pthread_barrier_wait(&barrier_division_static);

        //以线程数为步长，划分任务
        for (int i = k + threadID + 1; i < n; i += THREAD_NUM)
        {
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            a[i][k] = 0;
        }

        //做完消去，所有线程同步
        pthread_barrier_wait(&barrier_elimination_static);
    }
    pthread_exit(NULL);
}

//静态线程高斯消去
void gaussEliminationStatic(int n, float a[][MAX_N])
{
    //初始化barrier
    pthread_barrier_init(&barrier_division_static, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_elimination_static, NULL, THREAD_NUM);

    //创建线程
    pthread_t handle[THREAD_NUM];          //线程句柄
    threadParam_t threadParam[THREAD_NUM]; //线程参数
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_create(&handle[threadID], NULL, staticThreadFunc, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    //销毁barrier
    pthread_barrier_destroy(&barrier_division_static);
    pthread_barrier_destroy(&barrier_elimination_static);
}

//===按行划分的静态线程消去======================================================================================================================

//===按列划分（跳跃式）的静态线程消去======================================================================================================================

//定义barrier
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
        for (int j = division_start; j < n; j += THREAD_NUM)
        {
            a[k][j] /= a[k][k];
        }

        //所有线程同步
        pthread_barrier_wait(&barrier_division_static_onColumn_mode1);

        //以线程数为步长，划分任务
        a[k][k] = 1.0;
        for (int i = k + 1; i < n; i++)
        {
            int elimination_start = threadID + k + 1;
            for (int j = elimination_start; j < n; j += THREAD_NUM)
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
    //初始化barrier
    pthread_barrier_init(&barrier_division_static_onColumn_mode1, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_elimination_static_onColumn_mode1, NULL, THREAD_NUM);

    //创建线程
    pthread_t handle[THREAD_NUM];          //线程句柄
    threadParam_t threadParam[THREAD_NUM]; //线程参数
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_create(&handle[threadID], NULL, staticThreadFunc_onColumn_mode1, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    //销毁barrier
    pthread_barrier_destroy(&barrier_division_static_onColumn_mode1);
    pthread_barrier_destroy(&barrier_elimination_static_onColumn_mode1);
}

//===按列划分（跳跃式）的静态线程消去======================================================================================================================

//===按列划分（连续式）的静态线程消去======================================================================================================================

//定义barrier
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
        int extraTask_d = amount_d % THREAD_NUM;                                                       //额外的任务量
        int h_d = amount_d / THREAD_NUM;                                                               //步长
        int my_start_d = threadID * h_d + k + 1 + ((threadID < extraTask_d) ? threadID : extraTask_d); //当前线程的开始位置
        int my_end_d = my_start_d + h_d + ((threadID < extraTask_d) ? 1 : 0);                          //当前线程的结束位置
        // cout << "From Thead =" << threadID << " and k = " << k << " : " << my_start_d << " to " << my_end_d << endl;
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
            int extraTask_e = amount_e % THREAD_NUM;                                                       //额外的任务量
            int h_e = amount_e / THREAD_NUM;                                                               //步长
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
    //初始化barrier
    pthread_barrier_init(&barrier_division_static_onColumn_mode2, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_elimination_static_onColumn_mode2, NULL, THREAD_NUM);

    //创建线程
    pthread_t handle[THREAD_NUM];          //线程句柄
    threadParam_t threadParam[THREAD_NUM]; //线程参数
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_create(&handle[threadID], NULL, staticThreadFunc_onColumn_mode2, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    //销毁barrier
    pthread_barrier_destroy(&barrier_division_static_onColumn_mode2);
    pthread_barrier_destroy(&barrier_elimination_static_onColumn_mode2);
}

//===按列划分（连续式）的静态线程消去======================================================================================================================

//===SIMD并行化高斯消去算法======================================================================================================================

// SIMD化的高斯消去
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

//===SIMD并行化高斯消去算法======================================================================================================================

//===SIMD、Pthread并行化高斯消去算法======================================================================================================================

//定义barrier
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
        for (int i = k + 1 + threadID; i < n; i += THREAD_NUM)
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

// SIMD化的高斯消去
void gaussEliminationSIMD_Pthread(int n, float a[][MAX_N])
{
    //初始化barrier
    pthread_barrier_init(&barrier_division_static_SIMD, NULL, THREAD_NUM);
    pthread_barrier_init(&barrier_elimination_static_SIMD, NULL, THREAD_NUM);

    //创建线程
    pthread_t handle[THREAD_NUM];          //线程句柄
    threadParam_t threadParam[THREAD_NUM]; //线程参数
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        threadParam[threadID].threadID = threadID;
        threadParam[threadID].n = n;
        threadParam[threadID].a = (void *)a;
    }
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_create(&handle[threadID], NULL, staticThreadFunc_SIMD, &threadParam[threadID]);
    }

    //等待线程结束
    for (int threadID = 0; threadID < THREAD_NUM; threadID++)
    {
        pthread_join(handle[threadID], NULL);
    }

    //销毁barrier
    pthread_barrier_destroy(&barrier_division_static_SIMD);
    pthread_barrier_destroy(&barrier_elimination_static_SIMD);
}

//===SIMD、Pthread并行化高斯消去算法======================================================================================================================


//===计时函数======================================================================================================================
void getTime(int n, float a[][MAX_N], int mode, double &duration1)
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
    default:
        break;
    }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    duration1 += time_span.count();
    // cout<<"Time = "<<time_span.count()<<endl;
}

//===计时函数======================================================================================================================

//===主函数======================================================================================================================

int main()
{
    double duration1[7] = {0.0}; //持续时长
    int times = 20;              //重复测20次
    // int n = 512;                 //系数矩阵规模
    int size[8] = {4, 8, 16, 32, 64, 128, 256, 512};
    for (int i = 0; i < 8; i++)
    {
        int n = size[i];
        float a[n][MAX_N] = {0};
        float b[n][MAX_N] = {0};
        float c[n][MAX_N] = {0};
        float d[n][MAX_N] = {0};
        float e[n][MAX_N] = {0};
        float f[n][MAX_N] = {0};
        float g[n][MAX_N] = {0};

        int cnt = times;
        while (cnt--)
        {
            m_reset(n, a);
            matrixDeepCopy(n, b, a);
            matrixDeepCopy(n, c, a);
            matrixDeepCopy(n, d, a);
            matrixDeepCopy(n, e, a);
            matrixDeepCopy(n, f, a);
            matrixDeepCopy(n, g, a);

            getTime(n, a, SERIAL_FUNC, duration1[0]);
            getTime(n, b, DYNAMIC_FUNC, duration1[1]);
            getTime(n, c, STATIC_FUNC, duration1[2]);
            getTime(n, d, STATIC_FUNC_COLUMN_MODE1, duration1[3]);
            getTime(n, e, STATIC_FUNC_COLUMN_MODE2, duration1[4]);
            getTime(n, f, SIMD_FUNC, duration1[5]);
            getTime(n, g, STATIC_FUNC_SIMD, duration1[6]);
        }

        cout << endl
             << "n =           " << n << endl;
        cout << "Serial:       " << duration1[0] / times * 1000 << endl;
        cout << "Dynamic:      " << duration1[1] / times * 1000 << endl;
        cout << "Static:       " << duration1[2] / times * 1000 << endl;
        cout << "Static mode1: " << duration1[3] / times * 1000 << endl;
        cout << "Static mode2: " << duration1[4] / times * 1000 << endl;
        cout << "SIMD:         " << duration1[5] / times * 1000 << endl;
        cout << "SIMD_Pthread: " << duration1[6] / times * 1000 << "      end"endl;
    }

    // gaussEliminationSerial(n, a);
    // gaussEliminationStatic(n, b);
    // gaussEliminationDynamic(n, c);
    // gaussEliminationStatic_onColumn_mode1(n, d);
    // gaussEliminationStatic_onColumn_mode2(n, e);

    // cout << "Serial Gauss Elimination" << endl;
    // printMatrix(n, a);

    // cout << "Static Gauss Elimination" << endl;
    // printMatrix(n, b);

    // cout << "Dynamic Gauss Elimination" << endl;
    // printMatrix(n, c);

    // cout << "Static Gauss Elimination on Column mode1" << endl;
    // printMatrix(n, d);

    // cout << "Static Gauss Elimination on Column mode2" << endl;
    // printMatrix(n, e);

    return 0;
}