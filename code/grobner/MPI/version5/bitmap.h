/**
 * @file bitmap.h
 * @author GKC_NKCS (2012522@mail.nankai.edu.cn)
 * @brief
 * @version 0.2
 * @date 2022-06-26
 *
 * @copyright Copyright (c) 2022
 * @details this implements functions to do with bitmap and MPI mode
 *
 */
#include <iostream>
#include <string>
#include <sstream>
#include <omp.h>
using namespace std;
#define INDEX_BLOCK_SIZE 4
#define WORD_BITS 32
#define UNORDERED 0
#define ORDERED 1
#define NUM_THREADS 2
typedef struct BitManager
{
    int lftCol; // leftest column number
    int wrdLen; // words actually used
    int idxLen; // index length
    int *idx;   // index
    BitManager() : lftCol(-1), wrdLen(0), idxLen(0), idx(nullptr) {}
} BitManager;


void createBitMap(string *sparseLine, int *bitmap)
{
    if (*sparseLine == "")
        return;
    int value;
    int wrdIdx; // serial number of word
    int bitIdx; // serial number of bit
    stringstream ss(*sparseLine);
    while (ss >> value)
    {
        wrdIdx = value / WORD_BITS;
        bitIdx = value % WORD_BITS;
        *(bitmap + wrdIdx) |= (1 << bitIdx);
        // cout<<"value = "<<value<<" wrdIdx = "<<wrdIdx<<" bitIdx = "<<bitIdx<<endl;
    }
}

void createWnd(string *sparseWnd, int *wnd, int rows, int wrdLen, bool mode)
{
    if (mode == UNORDERED)
    {
        for (int i = 0; i < rows; i++)
        {
            // cout<<"sparse line: "<<*(sparseWnd+i)<<endl;
            // cout<<"i: "<<i<< endl;
            createBitMap(sparseWnd + i, wnd + wrdLen * i);
            // cout<<"I: "<<i<<" created"<<endl;
        }
        // cout<<"done"<<endl;
        return;
    }
    else
    {
        for (int i = 0; i < rows; i++)
        {
            stringstream ss(*(sparseWnd + i));
            int lftCol = -1;
            ss >> lftCol;
            // cout<<"sparse line: "<<*(sparseWnd+i)<<endl;
            // cout<<"lc: "<<lftCol<<" wrdLen: "<<wrdLen<<endl;
            createBitMap(sparseWnd + i, wnd + wrdLen * lftCol);
        }
    }
}

string toString(int *bitmap, int wrdLen)
{
    string result = "";
    stringstream ss;
    int value;
    int wrdIdx;
    int bitIdx;
    int flag = 0b10000000000000000000000000000000; // 1<<31
    for (int i = wrdLen - 1; i >= 0; i--)
    {
        if (*(bitmap + i) != 0)
        {
            wrdIdx = i;
            int tmp = *(bitmap + i);
            for (bitIdx = WORD_BITS - 1; bitIdx >= 0; bitIdx--, tmp <<= 1)
            {
                if ((tmp & flag) == 0)
                    continue;
                value = wrdIdx * WORD_BITS + bitIdx;
                ss << value << " ";
            }
        }
    }
    result.append(ss.str());
    return result;
}

string toString(int *bitmap, BitManager *bitManager)
{
    if (bitManager->lftCol == -1)
        return "";

    string result = "";
    stringstream ss;
    int lftCol;
    int wrdIdx;
    int bitIdx;
    int flag = 0b10000000000000000000000000000000;
    for (int i = bitManager->idxLen - 1; i >= 0; i--) // scan from tail to head
    {
        if (bitManager->idx[i] == 1) // check index (to accelerate)
        {
            for (int j = (i + 1) * INDEX_BLOCK_SIZE - 1; j >= i * INDEX_BLOCK_SIZE; j--)
            {
                if (*(bitmap + j) != 0) // check word
                {
                    wrdIdx = j;
                    int tmp = *(bitmap + j);
                    for (bitIdx = WORD_BITS - 1; bitIdx >= 0; bitIdx--, tmp <<= 1)
                    {
                        if ((tmp & flag) == 0)
                            continue;
                        lftCol = wrdIdx * WORD_BITS + bitIdx;
                        ss << lftCol << " ";
                    }
                }
            }
        }
    }
    result.append(ss.str());
    return result;
}

void toString(int *wnd, int wrdLen, string *result, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        result[i] = toString(wnd + i * wrdLen, wrdLen);
    }
}

void toString(int *wnd, int wrdLen, string *result, BitManager *bitManagers, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        result[i] = toString(wnd + i * wrdLen, bitManagers + i);
    }
}

void printWnd(int *wnd, int wrdLen, int rows)
{
    string *result = new string[rows];
    toString(wnd, wrdLen, result, rows);
    for (int i = 0; i < rows; i++)
    {
        cout << "Line " << i << " : " << result[i] << endl;
    }
    delete[] result;
    result = nullptr;
}

void printWnd(int *wnd, int wrdLen, BitManager *bitManagers, int rows)
{
    string *result = new string[rows];
    toString(wnd, wrdLen, result, bitManagers, rows);
    for (int i = 0; i < rows; i++)
    {
        cout << "Line " << i << " : " << result[i] << endl;
    }
    delete[] result;
    result = nullptr;
}

void xorBitmap(int *bitmap1, int *bitmap2, int wrdLen)
{
    for (int i = 0; i < wrdLen; i++)
    {
        *(bitmap1 + i) ^= *(bitmap2 + i);
    }
}

void buildBitManager(int *bitmap, int wrdLen, BitManager *bitManager)
{
    if (bitManager->idx == nullptr)
    {
        // cout << "null" << endl;
                    // cout << "wrdLen: " << wrdLen << endl;

        bitManager->idx = new int[wrdLen % 4 == 0 ? wrdLen / INDEX_BLOCK_SIZE : wrdLen / INDEX_BLOCK_SIZE + 1]{0};
    }else
    {
        for (int i = 0; i < (wrdLen % 4 == 0 ? wrdLen / 4 : wrdLen / 4 + 1); i++)
        {
            // cout << "wrdLen: " << wrdLen << endl;
            bitManager->idx[i] = 0;
            // cout << "i: " << i << endl;
        }
    }
    bitManager->idxLen = 0;
    bitManager->wrdLen = 0;
    bitManager->lftCol = -1;

    // cout<<toString(bitmap,wrdLen)<<endl;

    int flag = 0b10000000000000000000000000000000;
    for (int wrdIdx = wrdLen-1; wrdIdx >= 0; wrdIdx--) // scan from tail to head
    {
        if (*(bitmap + wrdIdx) != 0)
        {
            int tmp = *(bitmap + wrdIdx);
            for (int bitIdx = WORD_BITS - 1; bitIdx >= 0; bitIdx--, tmp <<= 1)
            {
                if ((tmp & flag) == 0)  continue;
                // cout<<"wrdIdx: "<<wrdIdx<<" BitIdx: "<<bitIdx<<endl;
                bitManager->lftCol = wrdIdx * WORD_BITS + bitIdx;
                bitManager->wrdLen = wrdIdx + 1;
                bitManager->idxLen = wrdIdx / INDEX_BLOCK_SIZE + 1;
                for (int i = 0; i < bitManager->idxLen; i++)
                {
                    for (int j = 0; j < INDEX_BLOCK_SIZE; j++)
                    {
                        if (*(bitmap + i * INDEX_BLOCK_SIZE + j) != 0)
                        {
                            bitManager->idx[i] = 1;
                            continue;
                        }
                    }
                }
                return;
            }
        }
    }
}

void buildBitManager(int *wnd, int wrdLen, BitManager *bitManager, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        // cout<<toString(wnd+wrdLen*i,wrdLen)<<endl;
        buildBitManager(wnd + i*wrdLen, wrdLen, bitManager + i);
        // cout<<"rows: "<<i<<endl;
        // cout<<toString(wnd+i,wrdLen)<<endl;
        // cout<<"lc: "<<bitManager[i].lftCol<<endl;
    }
}

void freeBitManager(BitManager *bitManager)
{
    delete[] bitManager->idx;
    bitManager->idx = nullptr;
    bitManager->idxLen = 0;
    bitManager->wrdLen = 0;
    bitManager->lftCol = -1;
}

void freeBitManager(BitManager *bitManagers, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        freeBitManager(bitManagers + i);
    }
}

void copyBitmapSingle(int *bitmap1, int *bitmap2, int wrdLen)
{
    for (int i = 0; i < wrdLen; i++)
    {
        *(bitmap2 + i) = *(bitmap1 + i);
    }
}

void copyBitMap(int *bitmap1, int *bitmap2, BitManager *bitManager1, BitManager *bitManager2)
{
    for (int i = 0; i < bitManager1->wrdLen; i++)
    {
        *(bitmap2 + i) = *(bitmap1 + i);
    }
    for (int i = 0; i < bitManager1->idxLen; i++)
    {
        bitManager2->idx[i] = bitManager1->idx[i];
    }
    bitManager2->idxLen = bitManager1->idxLen;
    bitManager2->wrdLen = bitManager1->wrdLen;
    bitManager2->lftCol = bitManager1->lftCol;
}

void xorBitmap(int *bitmap1, int *bitmap2, BitManager *bitManager1, BitManager *bitManager2)
{
    for (int i = 0; i < bitManager1->idxLen; i++)
    {
        if (bitManager1->idx[i] == 1 || bitManager2->idx[i] == 1) // check index
        {
            for (int j = i * INDEX_BLOCK_SIZE; j < (i + 1) * INDEX_BLOCK_SIZE; j++)
            {
                if (*(bitmap1 + j) != 0 || *(bitmap2 + j) != 0)
                {
                    *(bitmap1 + j) ^= *(bitmap2 + j);
                }
            }
        }
    }
    buildBitManager(bitmap1, bitManager1->wrdLen, bitManager1);
}