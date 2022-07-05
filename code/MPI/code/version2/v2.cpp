/**
 * @file v2.cpp
 * @author NKCS_GKC (2012522@mail.nankai.edu.cn)
 * @brief 1d row-looped distribution
 * @version 0.1
 * @date 2022-06-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include "mpi.h"
using namespace std;

#define n 1024 // for the convenience of programming, n is able to be divided by 16
#define block_size 64
#define PRECISION 999 // ensure precision of 0.0001

float M[n][n] = {0.0};
float M_B[n][n] = {0.0}; // backup of matrix

//===functions to do with matrix
void initMatrix(float a[][n]);               // init matrix
void copyMatrix(float a[][n], float b[][n]); // deep copy mtrix a to b
void printMatrix(float a[][n]);              // print matrix
void readMatrix(float a[][n]);               // read matrix from file saved
//===serial gauss elimination
void gaussEliminationSerial(float a[][n]);
//===gauss elimination implemented by MPI
void gaussEliminationMPI(float a[][n], int argc, char *argv[]);

int main(int argc, char *argv[])
{
    // initMatrix(M);
    // readMatrix(M);
    // printMatrix(M);
    // gaussEliminationSerial(M);
    gaussEliminationMPI(M, argc, argv);
    // printMatrix(M);
}

// init matirx
void initMatrix(float a[][n])
{
    srand((unsigned)time(0));

    // init matrix by random element in (-1,1)
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            float t = pow(-1, rand() % 2 + 1);
            float temp = rand() % (PRECISION) / (double)(PRECISION);
            a[i][j] = t * temp;
        }
    }

    // random linear combination
    for (int i = 0; i < n * 10; i++)
    {
        int i1 = rand() % n;
        int i2 = rand() % n;
        for (int j = 0; j < n; j++)
        {
            a[i1][j] += a[i2][j];
        }
    }
}

// deep copy a to b
void copyMatrix(float a[][n], float b[][n])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            b[i][j] = a[i][j];
        }
    }
}

// print matrix
void printMatrix(float a[][n])
{
    cout << "Printing Matrix, Lines " << n << endl;
    for (int i = 0; i < n; i++)
    {
        cout << "Line " << i << " : ";
        for (int j = 0; j < n; j++)
        {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// serial gauss elimination
void gaussEliminationSerial(float a[][n])
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

// read matrix initialized from file taged by n
void readMatrix(float a[][n])
{
    string matrixPath = "./" + to_string(n) + ".txt";
    ifstream ifStream;
    ifStream.open(matrixPath);
    if (!ifStream)
    {
        ifStream.close();
        initMatrix(a);
        ofstream ofStream(matrixPath);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                ofStream << a[i][j] << ' ';
                // cout<<a[i][j]<<" ";
            }
            ofStream << endl;
        }
        ofStream.close();
        ifStream.open(matrixPath);
    }
    int r = 0;
    stringstream stringStream;
    while (!ifStream.eof())
    {
        string line;
        getline(ifStream, line);
        stringStream << line;
        for (int i = 0; i < n; i++)
        {
            stringStream >> a[r][i];
        }
        r++;
    }
    ifStream.close();
}

void gaussEliminationMPI(float a[][n], int argc, char *argv[])
{
    int myid;              // rank of current processor
    int num;               // numbers of processors
    double s_time, e_time; // count time

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num);
    MPI_Status status;
    MPI_Datatype V;

    MPI_Type_vector(n / block_size / num, block_size * n, num * block_size * n, MPI_FLOAT, &V);
    MPI_Type_commit(&V);

    // if rank = 0, init matrix and distribute task, else receive data from processor 0
    if (myid == 0)
    {
        readMatrix(a);
        s_time = MPI_Wtime();
        for (int i = 1; i < num; i++)
        {
            int i_start = i * block_size;
            MPI_Send(&a[i_start][0], 1, V, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&a[myid * block_size][0], 1, V, 0, 0, MPI_COMM_WORLD, &status);
        // printMatrix(a);
    }

    // wait until each processor receive its data
    MPI_Barrier(MPI_COMM_WORLD);

    // for (int i = 0; i < num; i++)
    // {
    //     if (myid == i)
    //     {
    //     cout<<"myid = "<<myid<<endl;
    //         // cout << "myid: " << myid << endl;
    //         printMatrix(a);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    for (int k = 0; k < n; k++)
    {
        // corresponding processor does the division work and asks higher-ranked processors to do the elimination work together
        int root = k / block_size % num; // find the processor having Kth row
        if (myid == root)
        {
            for (int j = k + 1; j < n; j++)
            {
                a[k][j] /= a[k][k];
            }
            a[k][k] = 1;
        }
        MPI_Bcast(&a[k][0], n, MPI_FLOAT, root, MPI_COMM_WORLD);
        // if(myid == 0)
        // {
        //     cout<<"k = "<<k<<" root= "<<root<<endl;
        // }

        // [r_start, r_end)
        int r_start, r_end;
        if (myid == root)
        {
            r_start = k + 1;
            r_end = r_start % block_size == 0 ? r_start + block_size : r_start + block_size - r_start % block_size;
            if (myid != r_start / block_size % num)
            {
                r_start = r_start + (num - 1) * block_size;
                r_end = min(r_start + block_size, n);
            }
            while (r_start < n)
            {
                // cout << "k = " << k << " myid = " << myid << " r_start =  " << r_start << " r_end =  " << r_end << endl;
                for (int i = r_start; i < r_end; i++)
                {
                    for (int j = k + 1; j < n; j++)
                    {
                        a[i][j] -= a[i][k] * a[k][j];
                    }
                    a[i][k] = 0;
                }
                r_start = r_end - block_size + num * block_size;
                r_end += block_size * num;
            }
        }
        else
        {
            r_start = (myid < root ? (myid - root + num) : (myid - root)) % num * block_size + k / block_size * block_size;
            r_end = min(r_start + block_size, n);

            while (r_start < n)
            {
                // cout << "k = " << k << " myid = " << myid << " r_start =  " << r_start << " r_end =  " << r_end << endl;
                for (int i = r_start; i < r_end; i++)
                {
                    for (int j = k + 1; j < n; j++)
                    {
                        a[i][j] -= a[i][k] * a[k][j];
                    }
                    a[i][k] = 0;
                }
                r_start += num * block_size;
                r_end += block_size * num;
            }
        }
    }
    //     int root_1 = (k + 1) / block_size % num;
    //     int r_start, r_end;
    //     if (myid == root_1)
    //     {
    //         r_start = k + 1;
    //         r_end = ((k + 1) / block_size + 1) * block_size - 1;
    //     }
    //     else
    //     {
    //         r_start = ((k + 1) / block_size + (myid - root_1) % 4) * block_size;
    //         r_end = r_start + num * block_size;
    //     }
    //     while (r_start < n)
    //     {
    //         for (int i = r_start; i <= r_end; i++)
    //         {
    //             for (int j = k + 1; j < n; j++)
    //             {
    //                 a[i][j] -= a[i][k] * a[k][j];
    //             }
    //             a[i][k] = 0;
    //         }
    //         r_start = r_end;
    //         r_end += num * block_size;
    //     }

    //     // wait until this rank is finished
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    // send back to processor 0
    if (myid != 0)
    {
        MPI_Send(&a[myid * block_size][0], 1, V, 0, myid, MPI_COMM_WORLD);
    }
    else
    {
        for (int i = 1; i < num; i++)
        {
            int i_start = i * block_size;
            MPI_Recv(&a[i_start][0], 1, V, i, i, MPI_COMM_WORLD, &status);
        }
    }
    // if (myid == 0)
    // {
    //     printMatrix(a);
    // }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0)
    {
        e_time = MPI_Wtime();
        cout << "time consuming: " << e_time - s_time << endl;
    }
    MPI_Finalize();
}
