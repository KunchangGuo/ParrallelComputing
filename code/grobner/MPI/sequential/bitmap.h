/**
 * @file bitmap.h
 * @author GKC_NKCS (2012522@mail.nankai.edu.cn)
 * @brief
 * @version 0.1
 * @date 2022-06-26
 *
 * @copyright Copyright (c) 2022
 * @details this implements functions to do with bitmap
 *
 */
#include <iostream>
#include <string>
#include <sstream>
using namespace std;
#define INDEX_BLOCK_SIZE 4
#define WORD_BITS 32

/**
 * @brief BitMap defination
 * @param array (int *)mainbody of bitmap
 * @param index (int *)index to accelerate elimination
 * @param lc (int)leftest column number
 * @param wl (int)length of array
 * @param il (int)length of index
 */
typedef struct tagBitMap
{
    int *array;
    int *index;
    int lc;
    int wl;
    int il;
    tagBitMap() : array(nullptr), index(nullptr), lc(-1), wl(0), il(0) {}
} BitMap;

void createBitMap(string *sparseLine, BitMap *bitmap);
void createWnd(string *sparseMatrix, BitMap *bitmap, int n, bool flag);
string toString(BitMap *bitmap);
void toString(BitMap *bitmap, string result[], int n);
void xorBitmap(BitMap *bitmap1, BitMap *bitmap2);
void freeBitmap(BitMap *bitmap);
void freeBitmap(BitMap *bitmap, int n);
void copyBitmap(BitMap *bitmap1, BitMap *bitmap2);

/**
 * @brief convert sparse line to bitmap
 * @param sparseLine-(string*)pointer to single line read from file
 * @param bitmap (BitMap*)convertion result
 */
void createBitMap(string *sparseLine, BitMap *bitmap)
{
    // if sparseLine is null, return
    if (*sparseLine == "")
        return;

    // init essential params of bitmap
    stringstream ss(*sparseLine);
    int lc;
    int wl;
    int il;
    ss >> lc;
    wl = lc / WORD_BITS + 1;
    wl += (wl % INDEX_BLOCK_SIZE) == 0 ? 0 : (INDEX_BLOCK_SIZE - wl % INDEX_BLOCK_SIZE);
    il = wl / INDEX_BLOCK_SIZE;
    bitmap->wl = wl;
    bitmap->lc = lc;
    bitmap->il = il;
    bitmap->array = new int[wl]{0};
    bitmap->index = new int[il]{0};

    // set corresponding bits to 1
    int wn; // serial number of word
    int bn; // serial number of bit
    int in; // serial number of index block
    do
    {
        wn = lc / WORD_BITS;
        bn = lc % WORD_BITS;
        in = wn / INDEX_BLOCK_SIZE;
        bitmap->array[wn] |= (1 << bn);
        bitmap->index[in] = 1;
    } while (ss >> lc);
}

/**
 * @brief convert sparse matrix to bitmap
 *
 * @param sparseLine (string*) arrays of sparse lines
 * @param bitmap (BitMap*) convert result (array)
 * @param n size of window
 * @param flag choose mode   false:unsequential  true:sequential
 */
void createWnd(string *sparseMatrix, BitMap *bitmap, int n, bool flag)
{
    if (flag == false)
    {
        for (int i = 0; i < n; i++)
        {
            createBitMap(sparseMatrix + i, bitmap + i);
        }
    }
    else
    {
        int lc = -1;
        for (int i = 0; i < n; i++)
        {
            stringstream ss(sparseMatrix[i]);
            ss << sparseMatrix[i];
            ss >> lc;
            createBitMap(sparseMatrix + i, bitmap + lc);
        }
    }
}

/**
 * @brief string form of bitmap
 * @param bitmap (BitMap*)
 * @return string
 */
string toString(BitMap *bitmap)
{
    string result = "";
    if (bitmap == nullptr)
        return result;
    stringstream ss;
    int lc;
    int wn;
    int bn;
    for (int il = bitmap->il - 1; il >= 0; il--) // scan from tail to head
    {
        if (bitmap->index[il] != 0) // check index (to accelerate)
        {
            for (int wl = (il + 1) * INDEX_BLOCK_SIZE - 1; wl >= il * INDEX_BLOCK_SIZE; wl--)
            {
                if (bitmap->array[wl] != 0) // check word
                {
                    wn = wl;
                    int tmp = bitmap->array[wn];
                    int flag = 0b10000000000000000000000000000000;
                    for (bn = WORD_BITS - 1; bn >= 0; bn--, tmp <<= 1)
                    {
                        if ((tmp & flag) == 0)
                        {
                            continue;
                        }
                        lc = wn * WORD_BITS + bn;
                        ss << lc << " ";
                    }
                }
            }
        }
    }
    result.append(ss.str());
    return result;
}

/**
 * @brief convert bitmap matrix to sparse matrix
 *
 * @param bitmap (BitMap*) bitmap matrix
 * @param result (string[]) convert result
 * @param n (int) window size
 * @return string
 */
void toString(BitMap *bitmap, string result[], int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = toString(bitmap + i);
    }
}

/**
 * @brief print part of bitmaps
 * 
 * @param bitmap (BitMap*)
 * @param n (int) size of bitmaps
 */
void printWnd(BitMap* bitmap, int n)
{
    string *result = new string[n];
    toString(bitmap,result,n);
    for(int i=0;i<n;i++)
    {
        cout<<result[i]<<endl;
    }
    delete []result;
    result = nullptr;
}

/**
 * @brief bitmap1 = bitmap1 xor bitmap2
 *
 * @param bitmap1 (BitMap*)
 * @param bitmap2 (BitMap*)
 */
void xorBitmap(BitMap *bitmap1, BitMap *bitmap2)
{
    // execute elimination
    for (int il = 0; il < bitmap1->il; il++)
    {
        if (bitmap1->index[il] != 0 || bitmap2->index[il] != 0)
        {
            for (int wl = il * INDEX_BLOCK_SIZE; wl < (il + 1) * INDEX_BLOCK_SIZE; wl++)
            {
                if (bitmap1->array[wl] != 0 || bitmap2->array[wl] != 0)
                {
                    bitmap1->array[wl] ^= bitmap2->array[wl];
                }
            }
        }
    }

    // rebuild index and leftest column of bitmap1
    int current_il = 0;
    for (int il = 0; il < bitmap1->il; il++)
    {
        for (int wl = il * INDEX_BLOCK_SIZE; wl < (il + 1) * INDEX_BLOCK_SIZE; wl++)
        {
            if (bitmap1->array[wl] != 0)
            {
                bitmap1->index[il] = 1;
                current_il = il + 1;
                break;
            }
            bitmap1->index[il] = 0;
        }
    }
    bitmap1->il = current_il;
    bitmap1->wl = bitmap1->il * INDEX_BLOCK_SIZE;
    bitmap1->lc = -1;
    for (int il = bitmap1->il - 1; il >= 0; il--) // scan from tail to head
    {
        if (bitmap1->index[il] != 0) // check index (to accelerate)
        {
            for (int wl = (il + 1) * INDEX_BLOCK_SIZE - 1; wl >= il * INDEX_BLOCK_SIZE; wl--)
            {
                if (bitmap1->array[wl] != 0) // check word
                {
                    int tmp = bitmap1->array[wl];
                    int flag = 0b10000000000000000000000000000000;
                    for (int bn = WORD_BITS - 1; bn >= 0; bn--, tmp <<= 1)
                    {
                        if ((tmp & flag) == 0)
                        {
                            continue;
                        }
                        bitmap1->lc = wl * WORD_BITS + bn;
                        return;
                    }
                }
            }
        }
    }
}

/**
 * @brief release resource allocated by single existence bitmap
 *
 * @param bitmap (BitMap*)
 */
void freeBitmap(BitMap *bitmap)
{
    if (bitmap == nullptr)
        return;
    delete[] bitmap->array;
    delete[] bitmap->index;
    bitmap->array = nullptr;
    bitmap->index = nullptr;
    bitmap = nullptr;
}

/**
 * @brief release resource allocated by bitmaps
 *
 * @param bitmap (BitMap*[])
 * @param n size of bitmap
 */
void freeBitmap(BitMap *bitmap, int n)
{
    for (int i = 0; i < n; i++)
    {
        if ((bitmap + i) != nullptr)
            freeBitmap(bitmap + i);
    }
    delete[] bitmap;
    bitmap = nullptr;
}

/**
 * @brief create bitmap2 using bitmap1
 *
 * @param bitmap1 (BitMap*)
 * @param bitmap2 (BitMap*)
 */
void copyBitmap(BitMap *bitmap1, BitMap *bitmap2)
{
    string sparseLine = toString(bitmap1);
    createBitMap(&sparseLine, bitmap2);
}