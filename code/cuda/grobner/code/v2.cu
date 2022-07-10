/**
 * @file v2.cpp
 * @author GKC_NKCS (2012522@mail.nankai.edu.cn)
 * @brief
 * @version 0.2
 * @date 2022-07-09
 *
 * @copyright Copyright (c) 2022
 * @details mainbody of CUDA gauss elimination
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
int* table1;
int* table2;

int threadsPerBlock;
int blocksPerGrid;

string basePath = "F:/大二下课程/并行计算/期末研究报告相关材料/data/Groebner/";
//string basePath = "/home/bill/Desktop/para/src/Groebner/";
string examplePath = basePath + getExampleName(1);

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
    checkCuda(cudaFree(table1));
    checkCuda(cudaFree(table2));
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

    size_t table1Size = wndSize1 * sizeof(int);
    size_t table2Size = wndSize * sizeof(int);
    checkCuda(cudaMallocManaged(&table1, table1Size));
    checkCuda(cudaMallocManaged(&table2, table2Size));
    initWith<<<blocksPerGrid,threadsPerBlock>>>(table1, wndSize1, -1);
    initWith<<<blocksPerGrid,threadsPerBlock>>>(table2, wndSize, -1);
    checkCuda(cudaDeviceSynchronize());
}

__global__ 
void eliminationKernelFunction(int* eliminatant, int* eliminator, int* table1, int* table2, int rank, int wndSize1, int wrdLen)
{
    __shared__ int correspondingRow;
    int row = blockIdx.x + rank;
    int i = row;
    int j;
    while (i < wndSize1)
    {
        if(threadIdx.x == 0)
            correspondingRow = table1[i];
        __syncthreads();
        if (correspondingRow != -1)
        {
            for(j = threadIdx.x;j<wrdLen;j+=blockDim.x)
            {
                eliminatant[i * wrdLen + j] ^= eliminator[correspondingRow * wrdLen + j];
            }
        }
        if (threadIdx.x == 0)
        {
            table1[i] = -1;
            int flag = 0b10000000000000000000000000000000;
            bool out = true;
            for (int wrdIdx = wrdLen - 1; out && wrdIdx >= 0; wrdIdx--) // scan from tail to head
            {
                int tmp = *(eliminatant + i * wrdLen + wrdIdx);
                if (tmp != 0)
                {
                    for (int bitIdx = WORD_BITS - 1; bitIdx >= 0; bitIdx--, tmp <<= 1)
                    {
                        if ((tmp & flag) == 0)  continue;
                        table1[i] = wrdIdx * WORD_BITS + bitIdx;
                        out = false;
                        break;
                    }
                }
            }
            correspondingRow = table1[i];
        }
        __syncthreads();
        if (correspondingRow == -1 || table2[correspondingRow] == -1)
        {
            i += gridDim.x;
        }
    }
}

void gaussian()
{
    int i, j;
    int* tmp = new int[wrdLen] {0};
    BitManager tmpManager;
    for (i = 0; i < wndSize1; i++)
    {
        buildBitManager(tmp, wrdLen, &tmpManager);

        int tmpLftCol = tmpManager.lftCol;
        if (tmpLftCol != -1 && eliminatorManager[tmpLftCol].lftCol == -1)
        {
            copyBitMap(tmp, eliminator + tmpLftCol * wrdLen, &tmpManager, eliminatorManager + tmpLftCol);
        }

        for (int idx = 0; idx < wndSize; idx++)
        {
            table2[idx] = eliminatorManager[idx].lftCol;
        }

        for (j = i; j < wndSize1; j++)
        {
            int lftCol = eliminatantManager[j].lftCol;
            table1[j] = eliminatorManager[lftCol].lftCol == -1 ? -1 : lftCol;
        }
        eliminationKernelFunction<<<blocksPerGrid, threadsPerBlock, 64>>>(eliminatant, eliminator, table1, table2, i, wndSize1, wrdLen);
        checkCuda(cudaDeviceSynchronize());
        checkCuda(cudaGetLastError());

        buildBitManager(eliminatant + wrdLen * i, wrdLen, eliminatantManager + i, wndSize1 - i);
        copyBitmapSingle(eliminatant+i*wrdLen, tmp, wrdLen);
    }
}

void write()
{
    string* result = new string[wndSize1];
    toString(eliminatant, wrdLen, result, wndSize1);
    writeResult(examplePath, result, wndSize1);
}