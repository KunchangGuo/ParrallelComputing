#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <unordered_set>
#include <bitset>
#include <chrono>
#include <arm_neon.h>
#include <unistd.h>
#include <string.h>
#include <semaphore>
using namespace std;

#define THREAD_NUM 1

const int maxColunmAmount = 562;

class BitMap;
class GrobnerBasedGaussElimination;

typedef struct {
    int threadID;
    fstream* eliminatant;
    fstream* eliminator;
    vector<string>* eliminatantSparseWindow;
    vector<string>* eliminatorSparseWindow;
    GrobnerBasedGaussElimination* grobner; //类对象，这里传递this指针
}threadParam_t;

sem_t sem_main;
set_t sem_worker;

void *eliminateFunc(void* param);
void *xorFunc(void* param);

//===位图类===
class BitMap
{
private:
    int wordLength = 32; //表示一个“字”长，“字”组成一行位图，这里“字”使用int类型，所以默认是32

    bool *secondLevelIndex;         //二级索引
    int secondLevelIndexLength = 0; //表示二级索引的长度，即有多少个“字”

    int *bits;              //位图主体部分，记录各个位的数据
    int highestNumber = -1; //最左端有值的列号，可以用于快速判断是否处在消元范围内

    bool eliminatorFlag = false; //对于被消元子，可能存在与后读取的消元子重复消去的可能，因此需要特殊标记

public:
    BitMap(){};  //构造函数
    ~BitMap(){}; //析构函数

    //位图初始化，使用矩阵列数作为输入参数，初始化二级索引长度
    void init(int columnAmount = maxColunmAmount)
    {
        secondLevelIndexLength = columnAmount / wordLength + 1;   //根据列数计算需要多少个“字”
        secondLevelIndexLength += (4+THREAD_NUM) - secondLevelIndexLength % (4+THREAD_NUM); //调整数据长度为4的倍数
        bits = new int[secondLevelIndexLength]{0};                //初始化位图主体
        secondLevelIndex = new bool[secondLevelIndexLength]{0};   //将二级索引全部初始化为0，当索引块中值不为0时，将该索引快初始化为1；
    }

    //置位函数
    void setBit(int x)
    {
        int secondLevelPosition = x / wordLength;               //计算出置位点在第几个字处
        int firstLevelPosition = x % wordLength;                //计算出置位点在该字的第几个位上
        bits[secondLevelPosition] |= (1 << firstLevelPosition); //置位
    }

    //更新二级索引
    void refreshSecondLevelIndex()
    {
        //遍历每个块
        for (int i = 0; i < secondLevelIndexLength; i++)
        {
            //如果这个块有值，即大于0，则将对应二级索引置为1
            if (bits[i] != 0)
            {
                secondLevelIndex[i] = 1;
            }
        }
    }

    //更新最左端列号
    void refreshHighestNumber()
    {
        //首先遍历二级索引，当一个索引块有值时，进行深层次的遍历
        for (int i = secondLevelIndexLength - 1; i >= 0; i--)
        {
            //没有值，继续找下一个块儿
            if (secondLevelIndex[i] == 0)
            {
                continue;
            }
            //有值，继续往下找
            int tempElement = bits[i];                     // 暂存这一个块儿
            int flag = 0b10000000000000000000000000000000; //用以寻找第一个出现的列号
            //通过位运算找到有值的点
            for (int j = wordLength - 1; j >= 0; j--)
            {
                //第j位为0
                if ((tempElement & flag) == 0)
                {
                    tempElement <<= 1; //左移一位，以便下次检查右边的低下一位
                    continue;
                }
                //第j位不为0，更新后直接返回
                this->highestNumber = i * wordLength + j; //找到最左端的列号
                return;
            }
        }
        //如果找不到，则置为-1
        this->highestNumber = -1;
    }

    //标记当前行为升格后的消元子
    void setToEliminator()
    {
        this->eliminatorFlag = true;
    }

    //判断当前行是否是升格后的消元子
    bool isEliminator()
    {
        return this->eliminatorFlag;
    }

    //获取最左端列号
    int getHighestNumber()
    {
        return this->highestNumber;
    }

    //稀疏矩阵行向量到位图的转化
    void sparseRowToBitMap(string &sparseRowString)
    {
        //将字符串型行向量转化为数值型行向量
        stringstream stringStream(sparseRowString); //先读取稀疏行向量
        vector<int> tempSparseRow;                  //将字符串行向量转化为整数型行向量
        int tempRowIndex = 0;                       //读取时的临时索引
        int tempElement;                            //承接每个元素的临时变量
        while (stringStream >> tempElement)         //完成读入工作，转换到数值型的稀疏行向量
        {
            tempSparseRow.push_back(tempElement);
            tempRowIndex++;
        }

        //初始化bits主体，对稀疏向量的每个元素进行对应位图的置位
        for (int i = 0; i < tempSparseRow.size(); i++)
        {
            setBit(tempSparseRow[i]); //逐个遍历数值型稀疏矩阵完成置位
        }

        //更新二级索引
        refreshSecondLevelIndex();

        //更新最左端列号
        refreshHighestNumber();
    }

    //位图到稀疏矩阵行向量的转化，结果存储到targetSparseRow中
    void bitMapToSparseRow(vector<int> &targetSparseRow)
    {
        //首先遍历二级索引，当一个索引块有值时，进行深层次的遍历
        for (int i = 0; i < secondLevelIndexLength; i++)
        {
            //没有值，继续找下一个块儿
            if (secondLevelIndex[i] == 0)
            {
                continue;
            }
            //有值，继续往下找
            int tempElement = bits[i]; // 暂存这一个块儿
            //通过位运算找到有值的点
            for (int j = 0; j < wordLength; j++)
            {
                //第j位为0
                if ((tempElement & 1) == 0)
                {
                    tempElement >>= 1; //右移一位，以便下次检查左边的下一位
                    continue;
                }
                //第j位不为0
                targetSparseRow.push_back(i * wordLength + j);
                tempElement >>= 1;
            }
        }
    }

    //位图转字符串，用于输出
    string toString()
    {
        string resultString = "";
        vector<int> tempSparseRow;
        bitMapToSparseRow(tempSparseRow);
        stringstream stringStream;
        for (int i = tempSparseRow.size() - 1; i >= 0; i--)
        {
            stringStream << tempSparseRow[i] << " ";
        }
        resultString.append(stringStream.str());
        return resultString;
    }

    //两个位图的异或运算，这里是将消元子的行作为参数进行传递
    void xorNormal(BitMap &eliminator)
    {
        //遍历二级索引
        for (int i = 0; i < secondLevelIndexLength; i++)
        {
            //当指向的这个块不同时为0时，进行异或操作
            if (secondLevelIndex[i] != 0 || eliminator.secondLevelIndex[i] != 0)
            {
                bits[i] ^= eliminator.bits[i];
            }
        }

        //更新二级索引
        refreshSecondLevelIndex();

        //更新最左端列号
        refreshHighestNumber();
    }

    //SIMD化的异或运算
    void xorSIMD(BitMap &eliminator)
    {
        int32x4_t t1;
        int32x4_t t2;
        //对连续的所有块进行异或
        for (int i = 0; i < secondLevelIndexLength; i += 4)
        {
            t1 = vld1q_s32(bits + i);
            t2 = vld1q_s32(eliminator.bits + i);
            t1 = veorq_s32(t1, t2);
            vst1q_s32(bits + i, t1);
            vst1q_s32(eliminator.bits + i, t2);
        }

        //更新二级索引
        refreshSecondLevelIndex();

        //更新最左端列号
        refreshHighestNumber();
    }

    //按ID划分的异或运算
    void xorPthread(BitMap &eliminator, int threadID)
    {
        int myStart = secondLevelIndexLength / THREAD_NUM * threadID;
        int myEND = secondLevelIndexLength / THREAD_NUM * (threadID+1);
        for (int i = myStart; i < myEND; i++)
        {
            //当指向的这个块不同时为0时，进行异或操作
            if (secondLevelIndex[i] != 0 || eliminator.secondLevelIndex[i] != 0)
            {
                bits[i] ^= eliminator.bits[i];
            }
        }
    }

    //按ID划分的SIMD异或运算
    void xorSIMDandPthread(BitMap &eliminator, int threadID)
    {
        int myStart = secondLevelIndexLength / THREAD_NUM * threadID;
        int myEND = secondLevelIndexLength / THREAD_NUM * (threadID+1);
        int32x4_t t1;
        int32x4_t t2;
        //对任务范围内连续的块进行异或
        for (int i = myStart; i < myEND; i += 4)
        {
            t1 = vld1q_s32(bits + i);
            t2 = vld1q_s32(eliminator.bits + i);
            t1 = veorq_s32(t1, t2);
            vst1q_s32(bits + i, t1);
            vst1q_s32(eliminator.bits + i, t2);
        }
    }
};

//======================================================================================
//======================================================================================

//===特殊高斯消元类===
class GrobnerBasedGaussElimination
{
public:
    string eliminatantPath = "";      //被消元子源文件的绝对路径
    string eliminatorOriginPath = ""; //消元子源文件的绝对路径
    string resultPath = "";           //结果文件的绝对路径
    string eliminatorPath = "";       //消元子文件的绝对路径

    int eliminatantWindowSize = 12;   //一轮做消去的被消元子的行数，即滑动窗在被消元子上的大小
    int eliminatorWindowSize = 12;    //一轮做消去的消元子的行数，即滑动窗在消元子上的大小
    vector<BitMap> eliminatantWindow; //被消元子滑动窗
    vector<BitMap> eliminatorWindow;  //消元子滑动窗
    
    unordered_set<int> eliminatorLeftHashSet; //消元子最左端首元素的哈希集合
    unordered_set<int> eliminatorWindowRange; //一批消元子最左端首元素的哈希集合

    BitMap* eliminatantLine; //当前进行消元的被消元子指针
    BitMap* eliminatorLine; //当前进行消元的消元子指针
    GrobnerBasedGaussElimination() {}
    ~GrobnerBasedGaussElimination() {}

    //输入例子的根目录，根据根目录完成初始化、消去计算的工作
    void init(string basicDirectory)
    {
        //初始化文件绝对路径
        this->eliminatantPath = basicDirectory + "/被消元行.txt";
        this->eliminatorOriginPath = basicDirectory + "/消元子.txt";
        this->eliminatorPath = basicDirectory + "/eliminator.txt";
        this->resultPath = basicDirectory + "/result.txt";

        // //初始化窗口大小
        // setWindowSize(basicDirectory);

        //拷贝消元子源文件中的数据到消元子文件中
        copyEliminator();

        //清空结果文件和临时文件
        clearResultFile();

        //初始化消元子哈希索引
        initEliminatorHashSet();

        //问题求解
        solve();
    }

    void setWindowSize(string directoryPath)
    {
        cout << directoryPath << endl;
        int i1 = directoryPath.find("矩阵列数") + 8;
        int i2 = directoryPath.find("，非零消元子");
        int len1 = i2 - i1;
        int i3 = directoryPath.find("非零消元子") + 10;
        int i4 = directoryPath.find("，被消元行");
        int len2 = i4 - i3;
        int i5 = directoryPath.find("被消元行") + 8;
        int i6 = directoryPath.length();
        int len3 = i6 - i5;
        string columns = directoryPath.substr(i1, len1);
        string eliminatorRows = directoryPath.substr(i3, len2);
        string eliminatantRows = directoryPath.substr(i5, len3);
        stringstream ss;
        ss << columns << " " << eliminatorRows << " " << eliminatantRows;
        int temp1, temp2, temp3;
        ss >> temp1 >> temp2 >> temp3;
        cout << temp1 << " " << temp2 << " " << temp3 << endl;
        this->eliminatorWindowSize = temp2 < 2048 ? temp2 : 2048;
        this->eliminatantWindowSize = temp3 < 2048 ? temp3 : 2048;
    }

    void copyEliminator()
    {
        ifstream eliminatorOrigin(eliminatorOriginPath, ios_base::in | ios_base::binary);
        ofstream eliminatorDestination(eliminatorPath, ios_base::trunc | ios_base::out | ios_base::binary);
        char *buf = new char[1024];
        int n;
        while (!eliminatorOrigin.eof())
        {
            eliminatorOrigin.read(buf, 1024);
            n = eliminatorOrigin.gcount();
            eliminatorDestination.write(buf, n);
        }
        eliminatorOrigin.close();
        eliminatorDestination.close();
    }

    void clearResultFile()
    {
        fstream resultFile;
        resultFile.open(resultPath, ios::out | ios::trunc);
        resultFile.close();
    }

    void initEliminatorHashSet()
    {
        fstream eliminator(eliminatorPath, ios::in);
        stringstream ss;
        string line;
        int highestNumber;
        while (!eliminator.eof())
        {
            //首先读取整行
            getline(eliminator, line);
            ss << line;
            //将首个元素转化为整型存到集合中
            ss >> highestNumber;
            eliminatorLeftHashSet.insert(highestNumber);
            //清空stringstream，以便下次读取
            ss.str("");
        }
    }

    //对全局的被消元子进行消去
    void solve()
    {
        //逐批读入被消元子，对每一批被消元子进行消去，并将消去结果写入结果文件
        fstream eliminatant;
        eliminatant.open(eliminatantPath, ios::in);
        vector<string> eliminatantSparseWindow; //被消元子字符串型滑动窗
        fstream eliminator;
        eliminator.open(eliminatorPath, ios::in);
        vector<string> eliminatorSparseWindow;    //消元子字符串型滑动窗
        unordered_set<int> eliminatorWindowRange; //当前消元子滑动窗的范围

        //创建线程
        pthread_barrier_init(&barrier, NULL, THREAD_NUM);
        pthread_t handle[THREAD_NUM];
        threadParam_t threadParam[THREAD_NUM];
        for (int threadID = 0; threadID < THREAD_NUM; threadID++)
        {
            threadParam[threadID].threadID = threadID;
            threadParam[threadID].eliminatant = &eliminatant;
            threadParam[threadID].eliminator = &eliminator;
            threadParam[threadID].eliminatantSparseWindow = &eliminatantSparseWindow;
            threadParam[threadID].eliminatorSparseWindow = &eliminatorSparseWindow;
            threadParam[threadID].grobner = this;
        }
        for (int threadID = 0; threadID < THREAD_NUM; threadID++)
        {
            pthread_create(&handle[threadID], NULL, eliminateFunc, &threadParam[threadID]);
        }

        //等待线程结束
        for (int threadID = 0; threadID < THREAD_NUM; threadID++)
        {
            pthread_join(handle[threadID], NULL);
        }
        //销毁barrier
        pthread_barrier_destroy(&barrier);

        eliminatant.close();
        eliminator.close();
    }

    //将运算完的被消元子行写入结果文件
    void writeEliminatantWindowToResultFile(BitMap eliminatantWindowLine)
    {
        //打开结果文件
        ofstream resultFile(resultPath, ios::out | ios::app);

        resultFile << eliminatantWindowLine.toString() << endl;

        //关闭结果文件
        resultFile.close();
    }

    //将消元子写入临时文件
    void writeEliminatorWindowToTempFile(vector<BitMap> &eliminatorWindow, int originSize)
    {
        //打开消元子
        ofstream tempFile(eliminatorPath, ios::app);

        //将消元子窗口中的每一行写入临时文件
        for (int i = originSize; i < eliminatorWindow.size(); i++)
        {
            tempFile << eliminatorWindow[i].toString() << endl;
        }

        //关闭临时文件
        tempFile.close();
    }
};

int main()
{
    string basePath = "/home/bill/Desktop/para/src/Groebner/";
    string exampleDirectory = "测试样例3 矩阵列数562，非零消元子170，被消元行53";
    string exampleDirectoryPath = basePath + exampleDirectory;

    GrobnerBasedGaussElimination g;
    g.init(exampleDirectoryPath);

    return 0;
}

void* eliminateFunc(void* param)
{
    //获取相应参数
    threadParam_t* p = (threadParam_t*)param;
    int threadID = p->threadID;
    fstream* eliminatant = p->eliminatant;
    fstream* eliminator = p->eliminator;
    vector<string>* eliminatantSparseWindow = p->eliminatantSparseWindow;
    vector<string>* eliminatorSparseWindow = p->eliminatorSparseWindow;
    GrobnerBasedGaussElimination* grobner = p->grobner;

    int originSize;

    //执行消去过程
    while (!eliminatant->eof())
    {
        //每次读入一批被消元子，直到达到滑动窗的大小上限或者读到文件末尾
        for (int i = 0; i < grobner->eliminatantWindowSize && !eliminatant->eof(); i++)
        {
            //读取每一行稀疏向量形式的被消元子
            string eliminatantSparseLine;
            getline((*eliminatant), eliminatantSparseLine);
            //将读取到的压入滑动窗中
            if (eliminatantSparseLine != "\0")
            {
                eliminatantSparseWindow->push_back(eliminatantSparseLine);
            }
        }

        //将被消元子滑动窗转化为位图形式
        for (int i = 0; i < eliminatantSparseWindow->size(); i++)
        {
            BitMap eliminatantBitMap;
            eliminatantBitMap.init();
            eliminatantBitMap.sparseRowToBitMap((*eliminatantSparseWindow)[i]);
            grobner->eliminatantWindow.push_back(eliminatantBitMap);
        }

        //将被消元子字符串型滑动窗清空，节省内存
        vector<string>().swap(*(eliminatantSparseWindow));

        //用每一批消元子对被消元子滑动窗进行消去
        for (int i = 0; i < grobner->eliminatantWindow.size(); i++)
        {
            //每次读入一批消元子，直到这一行被消元子变成空行或者成为升格的消元子
            while (true)
            {
                //读入一批消元子
                for (int j = 0; j < grobner->eliminatorWindowSize; j++)
                {
                    //读取每一行稀疏向量形式的消元子
                    string eliminatorSparseLine;
                    if (eliminator->eof())
                    {
                        eliminator->close(); //这里使用seekg的方法调整读入位置会读取空行，因此重新打开
                        eliminator->open(grobner->eliminatorPath, ios::in);
                    }
                    getline((*eliminator), eliminatorSparseLine);

                    //将读取到的压入消元子滑动窗中
                    if (eliminatorSparseLine != "\0")
                    {
                        eliminatorSparseWindow->push_back(eliminatorSparseLine);
                    }
                }

                //将字符串型消元子滑动窗转化为BitMap形式
                for (int j = 0; j < eliminatorSparseWindow->size(); j++)
                {
                    BitMap eliminatorBitMap;
                    eliminatorBitMap.init();
                    eliminatorBitMap.sparseRowToBitMap((*eliminatorSparseWindow)[j]);
                    grobner->eliminatorWindow.push_back(eliminatorBitMap);
                    grobner->eliminatorWindowRange.insert(eliminatorBitMap.getHighestNumber());
                }

                //记录当前消元子滑动窗的大小
                originSize = grobner->eliminatorWindow.size();

                //清空字符串型消元子滑动窗，节省内存
                vector<string>().swap((*eliminatorSparseWindow));

                //先判断当前行是否可以进行消元
                while (grobner->eliminatorWindowRange.find(grobner->eliminatantWindow[i].getHighestNumber()) != grobner->eliminatorWindowRange.end())
                {
                    for (int j = 0; j < grobner->eliminatorWindow.size(); j++)
                    {
                        if (grobner->eliminatorWindow[j].getHighestNumber() == grobner->eliminatantWindow[i].getHighestNumber())
                        {
                            //标记当前的消元子和被消元子
                            grobner->eliminatantLine = &grobner->eliminatantWindow[i];
                            grobner->eliminatorLine = &grobner->eliminatorWindow[j];

                            //唤醒其他线程一起做消去
                            sem_
                            // grobner->eliminatantWindow[i].xorNormal(grobner->eliminatorWindow[j]);
                            // grobner->eliminatantWindow[i].xorSIMD(grobner->eliminatorWindow[j]);
                            grobner->eliminatantWindow[i].xorPthread(grobner->eliminatorWindow[j],threadID);

                            grobner->eliminatantWindow[i].refreshSecondLevelIndex();
                            grobner->eliminatantWindow[i].refreshHighestNumber();
                            break;
                        }
                    }
                }

                //判断当前行是否可以升格
                if (grobner->eliminatantWindow[i].getHighestNumber() != -1 && grobner->eliminatorLeftHashSet.find(grobner->eliminatantWindow[i].getHighestNumber()) == grobner->eliminatorLeftHashSet.end())
                {
                    //可以升格，加入到消元子中，同时更新消元子索引，同时表示该行运算完毕，写入结果文件
                    grobner->eliminatantWindow[i].setToEliminator();
                    grobner->eliminatorWindow.push_back(grobner->eliminatantWindow[i]);
                    grobner->eliminatorLeftHashSet.insert(grobner->eliminatantWindow[i].getHighestNumber());
                    grobner->writeEliminatorWindowToTempFile(grobner->eliminatorWindow, originSize);
                    grobner->writeEliminatantWindowToResultFile(grobner->eliminatantWindow[i]);
                }

                //如果是空行，写入到结果文件中
                if (grobner->eliminatantWindow[i].getHighestNumber() == -1)
                {
                    grobner->writeEliminatantWindowToResultFile(grobner->eliminatantWindow[i]);
                }

                vector<BitMap>().swap(grobner->eliminatorWindow);
                grobner->eliminatorWindowRange.clear();

                //空行或者为升格后的消元子则直接对下一行进行操作
                if (grobner->eliminatantWindow[i].getHighestNumber() == -1 || grobner->eliminatantWindow[i].isEliminator())
                {
                    break;
                }
            }
        }
    }
    pthread_exit(NULL);
}

void* xorFunc(void* parm)
{
    threadParam_t* p = (threadParam_t*)parm;
    GrobnerBasedGaussElimination* grobner = p->grobner;
    int threadID = p->threadID;
}