#include <iostream>
#include <arm_neon.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace std;

const int maxN = 640; // ϵ����������ģ

// �������:a������b
void matrixDeepCopy(int n, float b[][maxN], float a[][maxN])
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            b[i][j]=a[i][j];
        }
    }
}

// ���и�˹��ȥ
void gaussEliminationOrdinary(int n, float a[][maxN])
{
    for(int k=0;k<n;k++)
    {
        for(int j=k+1;j<n;j++)
        {
            a[k][j]/=a[k][k];
        }
        a[k][k]=1.0;
        for(int i=k+1;i<n;i++)
        {
            for(int j=k+1;j<n;j++)
            {
                a[i][j]-=a[i][k]*a[k][j];
            }
            a[i][k]=0;
        }
    }
}

// SIMD���ĸ�˹��ȥ
void gaussEliminationOptimisedUnmatched(int n,float a[][maxN])
{
    float32x4_t vt; // ����a[k][k]
    float32x4_t va; // ����a[k][j]~a[k][j+3]
    float32x4_t vaik; // ����a[i][k]
    float32x4_t vakj; // ����a[k][j]~a[k][j+3]
    float32x4_t vaij; // ����a[i][j]~a[i][j+3]
    float32x4_t vx; // ����(a[k][j]~a[k][j+3])*a[i][k]��ֵ

    int k;
    int i;
    int j;
    int offset;
    for(k=0;k<n;k++)
    {
        vt=vld1q_dup_f32(a[k]+k); // ��a[k][k]�浽vt���ĸ�ͨ����
        for(j=k+1;j+4<=n;j+=4)
        {
            va=vld1q_f32(a[k]+j); // ��a[k][j]~a[k][j+3]�浽va���ĸ�ͨ����
            va=vdivq_f32(va,vt); // a[k][j] /= a[k][k]
            vst1q_f32(a[k]+j,va); // ��������д�ص��ڴ���
        }
        offset=(n-k-1)%4;
        for(j=n-offset;j<n;j++)
        {
            a[k][j]/=a[k][k]; // δ�����Ԫ��
        }
        a[k][k]=1.0;
        for(i=k+1;i<n;i++)
        {
            vaik=vld1q_f32(a[i]+k); // ��a[i][k]�浽vaik���ĸ�ͨ����
            for(j=k+1;j+4<=n;j+=4)
            {
                vakj=vld1q_f32(a[k]+j); // ����a[k][j]~a[k][j+3]
                vaij=vld1q_f32(a[i]+j); // ����a[i][j]~a[i][j+3]
                vx=vmulq_f32(vakj,vaik); // ����(a[k][j]~a[k][j+3])*a[i][k]
                vaij=vsubq_f32(vaij,vx); // ����(a[i][j]~a[i][j+3])-(a[k][j]~a[k][j+3])*a[i][k]
                vst1q_f32(a[i]+j,vaij); // ��a[i][j]~a[i][j+3]д�ص��ڴ���
            }
            for(j=n-offset;j<n;j++)
            {
                a[i][j]-=a[k][j]*a[i][k]; // δ�����Ԫ��
            }
            a[i][k]=0;
        }
    }
}

// �����SIMD���ĸ�˹��ȥ
void gaussEliminationOptimisedMatched(int n,float a[][maxN])
{
    float32x4_t vt; // ����a[k][k]
    float32x4_t va; // ����a[k][j]~a[k][j+3]
    float32x4_t vaik; // ����a[i][k]
    float32x4_t vakj; // ����a[k][j]~a[k][j+3]
    float32x4_t vaij; // ����a[i][j]~a[i][j+3]
    float32x4_t vx; // ����(a[k][j]~a[k][j+3])*a[i][k]��ֵ

    int k;
    int i;
    int j;
    int start;
    int leftover = n%4;
    int end = n-leftover;
    int handle;
    for(k=0;k<n;k++)
    {
        vt=vld1q_dup_f32(a[k]+k); // ��a[k][k]�浽vt���ĸ�ͨ����
        start = k+1;
        handle=start%4;
        // ��ͷ
        switch(handle)
        {
            case 0:
                break; // ��ʼλ���Ѿ�����
            case 1:
                a[k][start]/=a[k][k];
                a[k][start+1]/=a[k][k];
                a[k][start+2]/=a[k][k];
                start+=3;
                break;
            case 2:
                a[k][start]/=a[k][k];
                a[k][start+1]/=a[k][k];
                start+=2;
                break;
            case 3:
                a[k][start]/=a[k][k];
                start+=1;
                break;
        }
        for(j=start;j+4<=end;j+=4)
        {
            va=vld1q_f32(a[k]+j); // ��a[k][j]~a[k][j+3]�浽va���ĸ�ͨ����
            va=vdivq_f32(va,vt); // a[k][j] /= a[k][k]
            vst1q_f32(a[k]+j,va); // ��������д�ص��ڴ���
        }
        // ȥβ��
        for(j=end;j<n;j++)
        {
            a[k][j]/=a[k][k]; // β�Ͳ���
        }
        a[k][k]=1.0;
        for(i=k+1;i<n;i++)
        {
            vaik=vld1q_f32(a[i]+k); // ��a[i][k]�浽vaik���ĸ�ͨ����
            start = k+1;
            handle=start%4;
            // ��ͷ
            switch(handle)
            {
                case 0:
                    break; // ��ʼλ���Ѿ�����
                case 1:
                    a[i][start]-=a[i][k]*a[k][start];
                    a[i][start+1]-=a[i][k]*a[k][start+1];
                    a[i][start+2]-=a[i][k]*a[k][start+2];
                    start+=3;
                    break;
                case 2:
                    a[i][start]-=a[i][k]*a[k][start];
                    a[i][start+1]-=a[i][k]*a[k][start+1];
                    start+=2;
                    break;
                case 3:
                    a[i][start]-=a[i][k]*a[k][start];
                    start+=1;
                    break;
            }
            for(j=start;j+4<=end;j+=4)
            {
                vakj=vld1q_f32(a[k]+j); // ����a[k][j]~a[k][j+3]
                vaij=vld1q_f32(a[i]+j); // ����a[i][j]~a[i][j+3]
                vx=vmulq_f32(vakj,vaik); // ����(a[k][j]~a[k][j+3])*a[i][k]
                vaij=vsubq_f32(vaij,vx); // ����(a[i][j]~a[i][j+3])-(a[k][j]~a[k][j+3])*a[i][k]
                vst1q_f32(a[i]+j,vaij); // ��a[i][j]~a[i][j+3]д�ص��ڴ���
            }
            //  ȥβ
            for(j=end;j<n;j++)
            {
                a[i][j]-=a[k][j]*a[i][k]; // δ�����Ԫ��
            }
            a[i][k]=0;
        }
    }
}


void printMatrix(int n, float a[][maxN])
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            cout<<a[i][j]<<' ';
        }
        cout<<endl;
    }
    cout<<endl;
}

void m_reset(int n, float a[][maxN])
{
	srand((unsigned)time(0));
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<i;j++)
        {
            a[i][j]=0;
        }
        a[i][i]=1.0;
        for(int j=i+1;j<n;j++)
        {
		int t=pow(-1,rand()%2+1);
            a[i][j]=t*rand()/float(RAND_MAX);
        }
    }
    for(int k=0;k<n;k++)
    {
        for(int i=k+1;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                a[i][j]+=a[k][j];
            }
        }
    }
}

int main()
{
    int n; // ����Ϊn*n
    int times = 3; // �ظ���������Ϊ3 
    int amount = 16; // �����ľ����ģ��16��

    /*int sizeN[] = {2,3,4,5,
	    6,7,8,12,
	    16,20,24,28,
	    32,46,47,48,
	    49,64,96,128,
	    196,256,512,1024};*/
    int sizeN[]={5,72,80,88,96,104,112,120,125,126,127,128,129,130,131,132};
    for(int i=0;i<amount;i++)
    {
        n=sizeN[i];

        float a[n][maxN] = {0};
        float b[n][maxN] = {0};
        float c[n][maxN] = {0};
        m_reset(n,a); // ����һ��a
        matrixDeepCopy(n,b,a); // bΪa�ĸ���
        matrixDeepCopy(n,c,a); // bΪa�ĸ���
        // printMatrix(n,a);

        // Linux�¸߾��ȼ�ʱ,��ÿ����ģ���ظ�5���������Լ��ټ������
        time_t dsec_a = 0, dsec_b = 0, dsec_c = 0;
        long long dnsec_a = 0, dnsec_b = 0, dnsec_c = 0;

        for(int j=0;j<times;j++)
        {
            // ����δ�����SIMD���ĸ�˹��ȥʱ��
            struct timespec sts,ets;
            timespec_get(&sts,TIME_UTC); // �����Ŀ�ʼʱ��
            gaussEliminationOptimisedUnmatched(n,a); // ����
            timespec_get(&ets,TIME_UTC); // �����Ľ���ʱ��
            time_t dsec = ets.tv_sec-sts.tv_sec;
            long dnsec = ets.tv_nsec-sts.tv_nsec;
            if(dnsec<0)
            {
                dsec--;
                dnsec+=1000000000ll;
            }
            // printMatrix(n,a);
            // printf ("%lld.%09llds\n",dsec,dnsec);
            dsec_a += dsec; // ��������ʱ������ۼ�
            dnsec_a += dnsec;

            // ����δ���еĸ�˹��ȥʱ��
            timespec_get(&sts,TIME_UTC); // �����Ŀ�ʼʱ��
            gaussEliminationOrdinary(n,b); // ����
            timespec_get(&ets,TIME_UTC); // �����Ľ���ʱ��
            dsec = ets.tv_sec-sts.tv_sec;
            dnsec = ets.tv_nsec-sts.tv_nsec;
            if(dnsec<0)
            {
                dsec--;
                dnsec+=1000000000ll;
            }
            // printMatrix(n,a);
            // printf ("%lld.%09llds\n",dsec,dnsec);
            dsec_b += dsec; // ��������ʱ������ۼ�
            dnsec_b += dnsec;

	  // ���������SIMD��˹��ȥʱ��
            timespec_get(&sts,TIME_UTC); // �����Ŀ�ʼʱ��
            gaussEliminationOptimisedMatched(n,c); // ����
            timespec_get(&ets,TIME_UTC); // �����Ľ���ʱ��
            dsec = ets.tv_sec-sts.tv_sec;
            dnsec = ets.tv_nsec-sts.tv_nsec;
            if(dnsec<0)
            {
                dsec--;
                dnsec+=1000000000ll;
            }
            // printMatrix(n,a);
            // printf ("%lld.%09llds\n",dsec,dnsec);
            dsec_c += dsec; // ��������ʱ������ۼ�
            dnsec_c += dnsec;

            // ���³�ʼ��a�����b����
            m_reset(n,a);
            matrixDeepCopy(n,b,a);
	  matrixDeepCopy(n,c,a);
        }

        // ���ƽ����ʱ
        dsec_a/=times;
        dsec_b/=times;
        dnsec_a/=times;
        dnsec_b/=times;
        dsec_c/=times;
        dnsec_c/=times;
        
        printf("ordinary     n= %d:  ordinary   = %11d.%09llds\n",n,dsec_b,dnsec_b);
        printf("optimised  n= %d: optimised = %11d.%09llds\n",n,dsec_a,dnsec_a);
        printf("matched    n= %d: optimised = %11d.%09llds\n",n,dsec_c,dnsec_c);
    }
    return 0;
}
