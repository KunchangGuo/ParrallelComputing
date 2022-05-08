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
//声明条件变量
pthread_cond_t cond_dynamic;
//声明互斥量的属性
pthread_mutexattr_t mutex_dynamic_attr;
//定义动态线程消去的行
int dynamic_row = 1;
//定义线程的触发变量
int count = -1;
//定义当前消去的轮次
int rank1 = 0;

//动态线程消去函数
void *dynamicThreadFunc(void *param)
{
    threadParam_t *p = (threadParam_t *)param;
    int threadID = p->threadID;                            //获取线程ID
    int n = p->n;                                          //获取矩阵规模
    float(*a)[MAX_N] = static_cast<float(*)[MAX_N]>(p->a); //获取矩阵
    int i;                                                 //每个线程的任务行
    int k;                                                 //消去的轮次
    pthread_mutex_lock(&mutex_dynamic);                    //动态线程消去时，动态线程需要加锁
    while (count > 0)
    {
        pthread_cond_wait(&cond_dynamic, &mutex_dynamic); //等待条件变量
        cout<<"Thread "<<threadID<<" woken"<<endl;
        count--;
        cout<<"dynamic_row:"<<dynamic_row<<endl;
        cout<<"rank1:"<<rank1<<endl;
        if (dynamic_row >= n)
        {
            //表示当前执行所有消去完毕
            if (rank1 >= n)
            {
                pthread_mutex_unlock(&mutex_dynamic); //解锁
                pthread_exit(NULL);
                return NULL;
            }
            //否则轮数和下一次开始消去的行递增
            rank1++;
            dynamic_row = rank1 + 1;
        }
        else
        {
            dynamic_row++;
        }
        i = dynamic_row; //设定当前的任务行
        //执行消去
        for (int j = rank1 + 1; j < n; j++)
        {
            a[i][j] -= a[i][rank1] * a[rank1][j];
        }
        a[i][rank1] = 0;
        pthread_mutex_unlock(&mutex_dynamic); //解锁
    }
}

//动态线程高斯消去
void gaussEliminationDynamic(int n, float a[][MAX_N])
{
    //初始化互斥量并设定它的属性
    pthread_mutexattr_init(&mutex_dynamic_attr);
    pthread_mutexattr_setrobust(&mutex_dynamic_attr, PTHREAD_MUTEX_ROBUST);
    pthread_mutex_init(&mutex_dynamic, &mutex_dynamic_attr);

    //初始化条件变量
    pthread_cond_init(&cond_dynamic, NULL);

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

    while (true)
    {
        pthread_mutex_lock(&mutex_dynamic); //动态线程消去时，动态线程需要加锁
        if (rank1 >= n)
        {
            //回收所有创建的线程并返回
            for (int threadID = 0; threadID < THREAD_NUM; threadID++)
            {
                pthread_join(handle[threadID], NULL);
            }
            pthread_mutex_unlock(&mutex_dynamic);
            //销毁互斥量
            pthread_mutex_destroy(&mutex_dynamic);
            //销毁条件变量
            pthread_cond_destroy(&cond_dynamic);
            //销毁互斥量属性
            pthread_mutexattr_destroy(&mutex_dynamic_attr);
            return;
        }
        for (int j = rank1 + 1; j < n; j++)
        {
            a[rank1][j] /= a[rank1][rank1];
        }
        a[rank1][rank1] = 1.0;
        count++;
        pthread_cond_signal(&cond_dynamic);   //唤醒消去线程
        pthread_mutex_unlock(&mutex_dynamic); //解锁
        cout << "here" << endl;
        cout << "rank" << rank1 << endl;
        cout<<"count:"<<count<<endl;
    }
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
