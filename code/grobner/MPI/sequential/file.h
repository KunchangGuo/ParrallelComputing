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
#define ELIMINATANT 0
#define ELIMINATOR 1

void getParam(string filePath, int &wndSize, int &rows);
void getSparseMatrix(string filePath, string *sparseMatrix, int n, int file);
void writeResult(string filePath, string *sparseMatrix, int n);

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
    // cout << wndSize << ' ' << wndSize2 << ' ' << wndSize1 << endl;
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
void getSparseMatrix(string filePath, string *sparseMatrix, int n, int file)
{
    if (file == ELIMINATANT)
    {
        filePath += "/被消元行.txt";
    }
    else if (file == ELIMINATOR)
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
    filePath += "/resultFile.txt";
    fstream fStream(filePath, ios::out | ios::trunc);
    for (int i = 0; i < n; i++)
    {
        fStream << sparseMatrix[i] << endl;
    }
    fStream.close();
}