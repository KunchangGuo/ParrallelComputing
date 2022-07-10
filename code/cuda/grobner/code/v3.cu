/**
 * @file v3.cpp
 * @author GKC_NKCS (2012522@mail.nankai.edu.cn)
 * @brief
 * @version 0.3
 * @date 2022-07-09
 *
 * @copyright Copyright (c) 2022
 * @details mainbody of MPI gauss elimination with streams
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

#define streamNumber 5

double s_time;    // start time
double e_time;    // end time
int* eliminatant; // eliminatant wnd
int* eliminator;  // eliminatant wnd
int* tmp;
int* table;
int* row;
int wndSize;      // max cols
int wndSize1;     // rows of eliminatant wnd
int wndSize2;     // rows of eliminator wnd
int wrdLen;       // cols per row
BitManager* eliminatantManager;
BitManager* eliminatorManager;

int threadsPerBlock;
int blocksPerGrid;

string basePath = "F:/大二下课程/并行计算/期末研究报告相关材料/data/Groebner/";
// string basePath = "/home/bill/Desktop/para/src/Groebner/";
string examplePath = basePath + getExampleName(5);

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

    for (int i = 0; i < wndSize1; i++)
    {
        delete[] eliminatantManager[i].idx;
        eliminatantManager[i].idx = nullptr;
    }
    for (int i = 0; i < wndSize; i++)
    {
        if (eliminatorManager[i].idx != nullptr)
        {
            delete[] eliminatorManager[i].idx;
            eliminatorManager[i].idx = nullptr;
        }
    }
    delete[] eliminatantManager;
    delete[] eliminatorManager;
    eliminatantManager = nullptr;
    eliminatorManager = nullptr;
    checkCuda(cudaFree(eliminatant));
    checkCuda(cudaFree(eliminator));
    /*checkCuda(cudaFree(tmp));
    checkCuda(cudaFree(table));
    checkCuda(cudaFree(row));*/
    return 0;
}

__global__ void initWith(int* wnd, int len, int num)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int grid = blockDim.x * gridDim.x;
    int i = tid;
    while (i < len)
    {
        wnd[i] = num;
        i += grid;
    }
}

void initWithCPU(int* wnd, int len, int num)
{
    for (int i = 0; i < len; i++)
    {
        wnd[i] = num;
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

    /*size_t tmpSize = wrdLen * sizeof(int);
    size_t tableSize = streamNumber * sizeof(int);
    checkCuda(cudaMallocManaged(&tmp, tmpSize));
    checkCuda(cudaMallocManaged(&table, tableSize));
    checkCuda(cudaMallocManaged(&row, tableSize));
    initWith << <blocksPerGrid, threadsPerBlock >> > (tmp, wrdLen, 0);
    initWith << <blocksPerGrid, threadsPerBlock >> > (table, streamNumber, -1);
    initWith << <blocksPerGrid, threadsPerBlock >> > (row, streamNumber, -1);
    checkCuda(cudaDeviceSynchronize());
    checkCuda(cudaGetLastError());*/
    tmp = new int[wrdLen];
    table = new int[streamNumber];
    row = new int[streamNumber];
    initWithCPU(tmp, wrdLen, 0);
    initWithCPU(table, streamNumber, -1);
    initWithCPU(row, streamNumber, -1);
}

__global__ void eliminationKernelFunction(int* bitmap1, int* bitmap2, int wrdLen)
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
        table[0] = lftCol;
        row[0] = i;
        while (lftCol != -1 && eliminatorManager[lftCol].lftCol != -1)
        {
            int cnt = 1;
            for (int j = i + 1; j < wndSize1 && cnt < streamNumber; j++)
            {
                int tmpLftCol = eliminatantManager[j].lftCol;
                if (tmpLftCol != -1 && eliminatantManager[tmpLftCol].lftCol != -1)
                {
                    table[cnt] = tmpLftCol;
                    row[cnt] = j;
                    //cout << "cnt: " << cnt << " j: "<<j<<" table[cnt]: " << table[cnt] << " row[cnt]: " << row[cnt]<<endl;
                    cnt++;
                }
            }
            for (int j = 0; j < cnt; ++j)
            {
                cudaStream_t stream;
                checkCuda(cudaStreamCreate(&stream));
                //cout << "j = " << j << " row[j] " << row[j] << " table[j] " << table[j] << endl;
                eliminationKernelFunction << <blocksPerGrid, threadsPerBlock >> > (eliminatant + row[j] * wrdLen, eliminator + table[j] * wrdLen, wrdLen);
                checkCuda(cudaStreamDestroy(stream));
            }
            checkCuda(cudaDeviceSynchronize());
            checkCuda(cudaGetLastError());
            for (int j = 0; j < cnt; j++)
            {
                buildBitManager(eliminatant + row[j] * wrdLen, wrdLen, eliminatantManager + row[j]);
            }
            lftCol = eliminatantManager[i].lftCol;
            table[0] = lftCol;
        }

        // store the result of current row
        copyBitmapSingle(eliminatant + i * wrdLen, tmp, wrdLen);
    }
}

void write()
{
    string* result = new string[wndSize1];
    toString(eliminatant, wrdLen, result, wndSize1);
    writeResult(examplePath, result, wndSize1);
    delete[] result;
    result = nullptr;
}