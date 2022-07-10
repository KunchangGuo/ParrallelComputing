/**
 * @file v1.cpp
 * @author GKC_NKCS (2012522@mail.nankai.edu.cn)
 * @brief
 * @version 0.1
 * @date 2022-07-06
 *
 * @copyright Copyright (c) 2022
 * @details mainbody of MPI gauss elimination
 *
 */
#include <stdio.h>
#include <string>
#include <iostream>
#include "file.h"
#include "bitmap.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "CheckCuda.cuh"

double s_time;    // start time
double e_time;    // end time
int* eliminatant; // eliminatant wnd
int* eliminator;  // eliminatant wnd
int wndSize;      // max cols
int wndSize1;     // rows of eliminatant wnd
int wndSize2;     // rows of eliminator wnd
int wrdLen;       // cols per row
BitManager* eliminatantManager;
BitManager* eliminatorManager;

int threadsPerBlock;
int blocksPerGrid;

string basePath = "F:/大二下课程/并行计算/期末研究报告相关材料/data/Groebner/";
//string basePath = "/home/bill/Desktop/para/src/Groebner/";
string examplePath = basePath + getExampleName(4);

void init();
void gaussian();
void write();

int main(int argc, char* argv[])
{
    float cudaElaspedTime = 0.0;
    cudaEvent_t cudaStart;
    cudaEvent_t cudaEnd;
    cudaEventCreate(&cudaStart);
    cudaEventCreate(&cudaEnd);

    getParam(examplePath, wndSize1, wndSize2, wndSize); // get size of wnd

    /* init wnd and relavant params */
    init();

    /* conduct elimination and timing*/
    cudaEventRecord(cudaStart, 0);
    gaussian();
    checkCuda(cudaEventRecord(cudaEnd, 0));
    checkCuda(cudaEventSynchronize(cudaEnd));
    checkCuda(cudaEventElapsedTime(&cudaElaspedTime, cudaStart, cudaEnd));
    printf("cuda elasped time: %f ms.\n", cudaElaspedTime);

    /*  gather and output result */
    write();

    delete[] eliminatantManager;
    delete[] eliminatorManager;
    eliminatantManager = nullptr;
    eliminatorManager = nullptr;
    checkCuda(cudaFree(eliminatant));
    checkCuda(cudaFree(eliminator));
    return 0;
}

__global__ void initWith(int* wnd, int len, int num)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int grid = blockDim.x * gridDim.x;
    int i = tid;
    while (i<len)
    {
        wnd[i] = num;
        i += grid;
    }
}

void init()
{
    wrdLen = wndSize / WORD_BITS + 1;

    size_t eliminatantSize = wndSize1 * wrdLen * sizeof(int);
    size_t eliminatorSize = wndSize * wrdLen * sizeof(int);
    checkCuda(cudaMallocManaged(&eliminatant, eliminatantSize));
    checkCuda(cudaMallocManaged(&eliminator, eliminatorSize));

    int deviceId;
    int numberOfSMs;
    cudaGetDevice(&deviceId);
    cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);
    threadsPerBlock = 1024;
    blocksPerGrid = 32 * numberOfSMs;

    initWith << <blocksPerGrid, threadsPerBlock >> > (eliminatant, wndSize1 * wrdLen, 0);
    initWith << <blocksPerGrid, threadsPerBlock >> > (eliminator, wndSize * wrdLen, 0);
    checkCuda(cudaDeviceSynchronize());

    string* eliminatantSparseWnd = new string[wndSize1];
    getSparseMatrix(examplePath, eliminatantSparseWnd, wndSize1, ELIMINATANT);
    createWnd(eliminatantSparseWnd, eliminatant, wndSize1, wrdLen, UNORDERED);
    delete[] eliminatantSparseWnd;
    eliminatantSparseWnd = nullptr;

    string* eliminatorSparseWnd = new string[wndSize2];
    getSparseMatrix(examplePath, eliminatorSparseWnd, wndSize2, ELIMINATOR);
    createWnd(eliminatorSparseWnd, eliminator, wndSize2, wrdLen, ORDERED);
    delete[] eliminatorSparseWnd;
    eliminatorSparseWnd = nullptr;

    eliminatantManager = new BitManager[wndSize1];
    eliminatorManager = new BitManager[wndSize];
    buildBitManager(eliminatant, wrdLen, eliminatantManager, wndSize1);
    buildBitManager(eliminator, wrdLen, eliminatorManager, wndSize);
}

__global__ 
void eliminationKernelFunction(int* bitmap1, int* bitmap2, int wrdLen)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int grid = blockDim.x * gridDim.x;
    int i = tid;
    while (i < wrdLen)
    {
        bitmap1[i] ^= bitmap2[i];
        i += grid;
    }
}

void gaussian()
{
    int* tmp = new int[wrdLen] {0};
    BitManager tmpManager;
    for (int i = 0; i < wndSize1; i++)
    {
        // dealing with result of last row
        buildBitManager(tmp, wrdLen, &tmpManager);
        int tmpLftCol = tmpManager.lftCol;
        if (tmpLftCol != -1 && eliminatorManager[tmpLftCol].lftCol == -1)
        {
            copyBitMap(tmp, eliminator + tmpLftCol * wrdLen, &tmpManager, eliminatorManager + tmpLftCol);
        }
        // dealing with current row
        int lftCol = eliminatantManager[i].lftCol;
        while (lftCol != -1 && eliminatorManager[lftCol].lftCol != -1)
        {
            eliminationKernelFunction << <blocksPerGrid, threadsPerBlock >> > (eliminatant + i * wrdLen, eliminator + lftCol * wrdLen, wrdLen);
            checkCuda(cudaDeviceSynchronize());
            checkCuda(cudaGetLastError());
            buildBitManager(eliminatant + i * wrdLen, wrdLen, eliminatantManager + i);
            lftCol = eliminatantManager[i].lftCol;
        }
        // store the result of current row
        copyBitmapSingle(eliminatant+i*wrdLen, tmp, wrdLen);
    }
}

void write()
{
    string* result = new string[wndSize1];
    toString(eliminatant, wrdLen, result, wndSize1);
    writeResult(examplePath, result, wndSize1);
}