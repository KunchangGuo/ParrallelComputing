#include <iostream>
//#include <arm_neon.h>
#include <stdlib.h>
#include <chrono>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
using namespace std;

#define THREAD_NUM 4 //线程数
#define MAX_N 640    //系数矩阵最大规模

void *dynamicThreadFunc(void *parm); //动态线程函数声明
void *staticThreadFunc(void *parm);  //静态线程函数声明

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
    //随机选择两行进行线性相加，重复2n次
    for (int i = 0; i < 2 * n; i++)
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

//线程数据结构定义
typedef struct
{
    int threadID; //线程ID
    int n;        //矩阵规模
    void *a;      //矩阵
} threadParam_t;

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
        cout<<"From Thread :"<<threadID<<endl;
        //做完除法，所有线程同步
        pthread_barrier_wait(&barrier_division_dynamic);
        cout<<"done first step!"<<endl;

        //所有线程从row中取值，并递增，需要用互斥量保护
        while (true)
        {
            //有保护地获取当前任务行
            pthread_mutex_lock(&mutex_dynamic);
            i = dynamic_row;
            cout<<"From Thread :"<<threadID<<" i = "<<i<<endl;
            dynamic_row++;
            if(i>=n) 
            {
                //下一次k取k+1，要从k+2开始消去
                dynamic_row = k+2;
                i = k+2;
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

int main()
{
    int n = 10; //系数矩阵规模
    float a[n][MAX_N] = {0};
    float b[n][MAX_N] = {0};
    float c[n][MAX_N] = {0};
    m_reset(n, a);
    matrixDeepCopy(n, b, a);
    matrixDeepCopy(n, c, a);

    cout << "a" << endl;
    printMatrix(n, a);

    cout << "b" << endl;
    printMatrix(n, b);

    cout << "c" << endl;
    printMatrix(n, c);

    gaussEliminationSerial(n, a);
    gaussEliminationStatic(n, b);
    gaussEliminationDynamic(n, c);

    cout << "Serial Gauss Elimination" << endl;
    printMatrix(n, a);

    cout << "Static Gauss Elimination" << endl;
    printMatrix(n, b);

    cout << "Dynamic Gauss Elimination" << endl;
    printMatrix(n, c);

    return 0;
}
