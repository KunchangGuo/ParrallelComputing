/**
 * @file gaussian.cpp
 * @author GKC_NKCS (2012522@mail.nankai.edu.cn)
 * @brief
 * @version 0.1
 * @date 2022-06-26
 *
 * @copyright Copyright (c) 2022
 * @details this implements gauss elimination based on grobner basis
 *
 */
#include <iostream>
#include <string>
#include "bitmap.h"
#include "file.h"
#include <chrono>
using namespace std;

string basePath = "F:/大二下课程/并行计算/期末研究报告相关材料/data/Groebner/";
// string basePath = "/home/bill/Desktop/para/src/Groebner/";
string examplePath = basePath + getExampleName(7);
int wndSize;         // max cols
int wndSize1;        // rows of eliminatant
int wndSize2;        // rows of eliminator
BitMap *eliminatant; // to be eliminated
BitMap *eliminator;  // to eliminate
int *table;          // table storing leftest column number

void initEliminatant();
void initEliminator();
void gaussian();

int main()
{
    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();
    // init
    getParam(basePath + examplePath, wndSize1, wndSize2, wndSize);

    high_resolution_clock::time_point i_start = high_resolution_clock::now();
    initEliminatant();
    initEliminator();
    high_resolution_clock::time_point i_end = high_resolution_clock::now();
    duration<double> i_time_span = duration_cast<duration<double>>(i_end - i_start);
    cout<<"time: " << i_time_span.count()<<endl;

    // calculate
    gaussian();

    // output
    string *sparseMatrix = new string[wndSize1];
    toString(eliminatant, sparseMatrix, wndSize1);
    writeResult(basePath + examplePath, sparseMatrix, wndSize1);

    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end - start);
    cout<<"time: " << time_span.count()<<endl;

    cout<<"ratio: "<<endl<<"init: "<<i_time_span.count()/time_span.count()<<endl<<"calculate: "<<1-i_time_span.count()/time_span.count()<<endl;

    // free resource
    delete[] table;
    delete[] sparseMatrix;
    table = nullptr;
    sparseMatrix = nullptr;
    freeBitmap(eliminatant, wndSize1);
    freeBitmap(eliminator, wndSize);
    return 0;
}

/**
 * @brief init eliminatant from file
 *
 */
void initEliminatant()
{
    string *sparseMatrix = new string[wndSize1];
    getSparseMatrix(basePath + examplePath, sparseMatrix, wndSize1, ELIMINATANT);
    eliminatant = new BitMap[wndSize1];
    createWnd(sparseMatrix, eliminatant, wndSize1, false);
    delete[] sparseMatrix;
    sparseMatrix = nullptr;
}

/**
 * @brief init eliminator from file and set table
 *
 */
void initEliminator()
{
    string *sparseMatrix = new string[wndSize2];
    getSparseMatrix(basePath + examplePath, sparseMatrix, wndSize2, ELIMINATOR);
    eliminator = new BitMap[wndSize];
    createWnd(sparseMatrix, eliminator, wndSize2, true);
    string *resultMatrix = new string[wndSize];
    table = new int[wndSize];
    for (int i = 0; i < wndSize; i++)
    {
        table[i] = eliminator[i].lc;
    }
    delete[] sparseMatrix;
    sparseMatrix = nullptr;
}

/**
 * @brief elimination mainbody
 *
 */
void gaussian()
{
    for (int i = 0; i < wndSize1; i++)
    {
        int lc = eliminatant[i].lc;
        while (lc != -1 && table[lc] != -1)
        {
            xorBitmap(eliminatant + i, eliminator + lc);
            lc = eliminatant[i].lc;
        }
        if (table[lc] == -1)
        {
            copyBitmap(eliminatant + i, eliminator + lc);
            table[lc] = lc;
            continue;
        }
    }
}