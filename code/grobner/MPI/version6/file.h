/**
 * @file file.h
 * @author GKC_NKCS (2012522@mail.nankai.edu.cn)
 * @brief
 * @version 0.1
 * @date 2022-06-26
 *
 * @copyright Copyright (c) 2022
 * @details this implements functions to do with IO
 *
 */
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
using namespace std;
#define ELIMINATANT false
#define ELIMINATOR true


/**
 * @brief Get the wndSize and rows adjustively
 *
 * @param filePath example directory
 * @param wndSize size of eliminatant
 * @param rows size of eliminator
 */
void getParam(string filePath, int &wndSize1, int &wndSize2, int &wndSize)
{
    string paramPath = filePath + "/param.txt";
    fstream param(paramPath, ios::in);
    param >> wndSize;
    param >> wndSize2;
    param >> wndSize1;
    // cout << "Wndsize: " << wndSize << " wndSize2: " << wndSize2 << " WndSize1: " << wndSize1 << endl;
    param.close();
}

/**
 * @brief get sparse matrix from file
 *
 * @param filePath example directory path
 * @param sparseMatrix result matrix
 * @param n size of wnd
 * @param file determine eliminatant or eliminator to be read
 */
void getSparseMatrix(string filePath, string *sparseMatrix, int n, int mode)
{
    if (mode == ELIMINATANT)
    {
        filePath += "/被消元行.txt";
    }
    else if (mode == ELIMINATOR)
    {
        filePath += "/消元子.txt";
    }
    fstream fStream(filePath, ios::in);
    if (!fStream.eof())
    {
        for (int i = 0; i < n; i++)
        {
            getline(fStream, sparseMatrix[i]);
        }
    }
    fStream.close();
}

/**
 * @brief write result to file
 *
 * @param filePath example directory path
 * @param sparseMatrix elimination result
 * @param n wnd size
 */
void writeResult(string filePath, string *sparseMatrix, int n)
{
    filePath += "/resultFile6.txt";
    fstream fStream(filePath, ios::out | ios::trunc);
    for (int i = 0; i < n; i++)
    {
        fStream << sparseMatrix[i] << endl;
    }
    fStream.close();
}

string getExampleName(int number)
{
    switch (number)
    {
    case 1:
        return "测试样例1 矩阵列数130，非零消元子22，被消元行8";
    case 2:
        return "测试样例2 矩阵列数254，非零消元子106，被消元行53";
    case 3:
        return "测试样例3 矩阵列数562，非零消元子170，被消元行53";
    case 4:
        return "测试样例4 矩阵列数1011，非零消元子539，被消元行263";
    case 5:
        return "测试样例5 矩阵列数2362，非零消元子1226，被消元行453";
    case 6:
        return "测试样例6 矩阵列数3799，非零消元子2759，被消元行1953";
    case 7:
        return "测试样例7 矩阵列数8399，非零消元子6375，被消元行4535";
    case 8:
        return "测试样例8 矩阵列数23045，非零消元子18748，被消元行14325";
    case 9:
        return "测试样例9 矩阵列数37960，非零消元子29304，被消元行14921";
    case 10:
        return "测试样例10 矩阵列数43577，非零消元子39477，被消元行54274";
    case 11:
        return "测试样例11 矩阵列数85401，非零消元子5724，被消元行756";
    default:
        return "";
    }
}