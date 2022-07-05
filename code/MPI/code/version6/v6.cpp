/**
 * @file v6.cpp
 * @author NKCS_GKC (2012522@mail.nankai.edu.cn)
 * @brief pipeline gauss elimination + openmp
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
#include <omp.h>
using namespace std;

#define n 2048
#define precision 999 // ensure precision of 0.0001
#define NUM_THREADS 2

float (*M)[n];
// float M_B[n][n] = {0.0}; // backup of matrix

//===functions to do with matrix
void initMatrix(float a[][n]);               // init matrix
void copyMatrix(float a[][n], float b[][n]); // deep copy mtrix a to b
void printMatrix(float a[][n]);              // print matrix
void readMatrix(float a[][n]);               // read matrix from file saved
//===serial gauss elimination
void gaussEliminationSerial(float a[][n]);
//===gauss elimination implemented by MPI
void gaussEliminationMPIPipelineOpenMP(float a[][n], int argc, char *argv[]);

int main(int argc, char *argv[])
{
    // initMatrix(M);
    // readMatrix(M);
    // printMatrix(M);
    // gaussEliminationSerial(M);
    gaussEliminationMPIPipelineOpenMP(M, argc, argv);
    // printMatrix(M);
}

// init matirx
void initMatrix(float (*a)[n])
{
    srand((unsigned)time(0));

    // init matrix by random element in (-1,1)
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            float t = pow(-1, rand() % 2 + 1);
            float temp = rand() % (precision) / (double)(precision);
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
void printMatrix(float (*a)[n])
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
void gaussEliminationSerial(float (*a)[n])
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
void readMatrix(float (*a)[n])
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

void gaussEliminationMPIPipelineOpenMP(float (*a)[n], int argc, char *argv[])
{
    int myid;              // rank of current processor
    int num;               // numbers of processors
    double s_time, e_time; // count time
    float(*sub)[n];        // submatrix of each processor in 1d
    float *tmp;            // get data of row k
    int np;                // submatrix size of each processor
    int provided;          // thread safety level provided

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_FUNNELED)
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num);
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Status status;

    a = new float[n][n];
    np = n / num;          // lines of submatrix
    sub = new float[n][n]; // store tmp result during computation
    tmp = new float[n];    // store row k during computation

    if (myid == 0)
    {
        readMatrix(a);
        s_time = MPI_Wtime();
    }

    // distribute task to each processor
    MPI_Scatter(&a[0][0], n * np, MPI_FLOAT, &sub[0][0], n * np, MPI_FLOAT, 0, com);

    int i, j, row;
#pragma omp parallel private(i, j, row)
    for (i = 0; i < myid * np; i++)
    {
#pragma omp single
        {
            // Get data from last level of pipeline
            MPI_Recv(tmp, n, MPI_FLOAT, myid - 1, 0, com, &status);
            // send data to next level of pipeline
            if (myid != (num - 1))
            {
                MPI_Send(tmp, n, MPI_FLOAT, myid + 1, 0, com);
            }
        }

#pragma omp for
        // eliminate individually
        for (row = 0; row < np; row++)
        {
            // currently, i equals k, standing for the rank of elimination
            for (j = i + 1; j < n; j++)
            {
                sub[row][j] -= sub[row][i] * tmp[j];
            }
            sub[row][i] = 0;
        }
    }

#pragma omp parallel private(row, i, j)
    for (row = 0; row < np; row++)
    {
#pragma omp for
        // calculation first
        for (j = myid * np + row + 1; j < n; j++)
        {
            sub[row][j] = sub[row][j] / sub[row][np * myid + row];
        }
        sub[row][np * myid + row] = 1;

#pragma omp for
        for (i = 0; i < n; i++)
        {
            tmp[i] = sub[row][i];
        }
#pragma omp single
        // send result to next pipeline
        if (myid != (num - 1))
        {
            MPI_Send(tmp, n, MPI_FLOAT, myid + 1, 0, com);
        }

// eliminate individually
#pragma omp for
        for (i = row + 1; i < np; i++)
        {
            for (j = myid * np + row + 1; j < n; j++)
            {
                sub[i][j] -= sub[i][myid * np + row] * tmp[j];
            }
            sub[i][myid * np + row] = 0;
        }
    }

    // Synchronize and gather rows
    MPI_Barrier(com);
    MPI_Gather(&sub[0][0], n * np, MPI_FLOAT, &a[0][0], n * np, MPI_FLOAT, 0, com);

    // // to correct result
    // if (myid == 0)
    // {
    //     printMatrix(a);
    // }

    if (myid == 0)
    {
        e_time = MPI_Wtime();
        cout << "time consuming: " << e_time - s_time << endl;
    }

    MPI_Finalize();
}
