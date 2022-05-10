#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <unordered_set>
#include <bitset>
#include <chrono>
#include <arm_neon.h>
using namespace std;

int maxColunmAmount = 562;

//===位图类===
class BitMap
{
private:
	int wordLength = 32; //表示一个“字”长，“字”组成一行位图，这里“字”使用int类型，所以默认是32

	bool* secondLevelIndex;  //二级索引
	int secondLevelIndexLength = 0; //表示二级索引的长度，即有多少个“字”

	int* bits;       //位图主体部分，记录各个位的数据
	int highestNumber = -1; //最左端有值的列号，可以用于快速判断是否处在消元范围内

	bool eliminatorFlag = false; //对于被消元子，可能存在与后读取的消元子重复消去的可能，因此需要特殊标记

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

public:
	BitMap() {};  //构造函数
	~BitMap() {}; //析构函数

	//位图初始化，使用矩阵列数作为输入参数，初始化二级索引长度
	void init(int columnAmount = maxColunmAmount)
	{
		secondLevelIndexLength = columnAmount / wordLength + 1; //根据列数计算需要多少个“字”
		secondLevelIndexLength += 4 - secondLevelIndexLength % 4; //调整数据长度为4的倍数
		bits = new int[secondLevelIndexLength] {0};              //初始化位图主体
		secondLevelIndex = new bool[secondLevelIndexLength] {0}; //将二级索引全部初始化为0，当索引块中值不为0时，将该索引快初始化为1；
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
	void sparseRowToBitMap(string& sparseRowString)
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

			////===测试代码
			//cout << "Element: " << tempElement << endl;
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
	void bitMapToSparseRow(vector<int>& targetSparseRow)
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
	void xorNormal(BitMap& eliminator)
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

	void xorSIMD(BitMap& eliminator)
	{
		int32x4_t t1;
		int32x4_t t2;
		//对连续的所有块进行异或
		for(int i=0;i<secondLevelIndexLength;i+=4)
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
};

//===特殊高斯消元类===
class GrobnerBasedGaussElimination
{
private:
	string eliminatantPath = "";      //被消元子源文件的绝对路径
	string eliminatorOriginPath = ""; //消元子源文件的绝对路径
	string resultPath = "";           //结果文件的绝对路径
	string eliminatorPath = "";       //消元子文件的绝对路径

	int eliminatantWindowSize = 2048;   //一轮做消去的被消元子的行数，即滑动窗在被消元子上的大小
	int eliminatorWindowSize = 2048;    //一轮做消去的消元子的行数，即滑动窗在消元子上的大小
	int columns = 0;
	vector<BitMap> eliminatantWindow; //被消元子滑动窗
	vector<BitMap> eliminatorWindow;  //消元子滑动窗

	//消元子最左端首元素的哈希集合
	unordered_set<int> eliminatorLeftHashSet;

public:
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

		//初始化滑动窗大小
		setWindowSize(basicDirectory);

		//拷贝消元子源文件中的数据到消元子文件中
		copyEliminator();

		//清空结果文件和临时文件
		clearResultFile();

		//初始化消元子哈希索引
		initEliminatorHashSet();

		//问题求解
		solve();
	}

private:
	void setWindowSize(string directoryPath)
	{
		//根据形如"测试样例3 矩阵列数562，非零消元子170，被消元行53"生成窗口大小
		int i1 = directoryPath.find("矩阵列数") + 8;
		int i2 = directoryPath.find("，非零消元子");
		int len1 = i2 - i1;
		int i3 = directoryPath.find("非零消元子") + 10;
		int i4 = directoryPath.find("，被消元行");
		int len2 = i4 - i3;
		int i5 = directoryPath.find("被消元行") + 8;
		int i6 = directoryPath.length() ;
		int len3 = i6 - i5;
		string columns = directoryPath.substr(i1,len1);
		string eliminatorRows = directoryPath.substr(i3,len2);
		string eliminatantRows = directoryPath.substr(i5,len3);
		stringstream ss;
		ss << columns << " " << eliminatorRows << " " << eliminatantRows;
		int temp1, temp2, temp3;
		ss >> temp1 >> temp2 >> temp3;
		this->eliminatantWindowSize = temp3 < 2048 ? temp3 : 2048;
		this->eliminatorWindowSize = temp2 < 2048 ? temp2 : 2048;
		this->columns = temp1;
	}

	void copyEliminator()
	{
		ifstream eliminatorOrigin(eliminatorOriginPath, ios_base::in | ios_base::binary);
		ofstream eliminatorDestination(eliminatorPath, ios_base::trunc | ios_base::out | ios_base::binary);
		char* buf = new char[1024];
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
		while (true)
		{
			//每次读入一批被消元子，直到达到滑动窗的大小上限或者读到文件末尾
			for (int i = 0; i < eliminatantWindowSize && !eliminatant.eof(); i++)
			{
				//读取每一行稀疏向量形式的被消元子
				string eliminatantSparseLine;
				getline(eliminatant, eliminatantSparseLine);
				//将读取到的压入滑动窗中
				if (eliminatantSparseLine != "\0")
				{
					eliminatantSparseWindow.push_back(eliminatantSparseLine);
				}
			}

			//将被消元子滑动窗转化为位图形式
			vector<BitMap> eliminatantWindow; //被消元子滑动窗
			for (int i = 0; i < eliminatantSparseWindow.size(); i++)
			{
				BitMap eliminatantBitMap;
				eliminatantBitMap.init();
				eliminatantBitMap.sparseRowToBitMap(eliminatantSparseWindow[i]);
				eliminatantWindow.push_back(eliminatantBitMap);
			}

			//将被消元子字符串型滑动窗清空，节省内存
			vector<string>().swap(eliminatantSparseWindow);

			//将这被消元子滑动窗用所有消元子滑动窗过一遍
			eliminate(eliminatantWindow);

			if (eliminatant.eof())
			{
				eliminatant.close();
				return;
			}
		}
	}

	//对被消元子滑动窗进行消去
	void eliminate(vector<BitMap>& eliminatantWindow)
	{
		//逐批读入消元子，用每一批消元子对被消元子滑动窗进行消去
		fstream eliminator;
		eliminator.open(eliminatorPath, ios::in);
		vector<string> eliminatorSparseWindow; //消元子字符串型滑动窗
		vector<BitMap> eliminatorWindow; //消元子滑动窗
		unordered_set<int> eliminatorWindowRange; //当前消元子滑动窗的范围

		for (int i = 0; i < eliminatantWindow.size(); i++)
		{
			//每次读入一批消元子，直到这一行被消元子变成空行或者成为升格的消元子
			while (true)
			{
				//空行或者为升格后的消元子则直接进行下一步
				if (eliminatantWindow[i].getHighestNumber() == -1 || eliminatantWindow[i].isEliminator())
				{
					break;
				}

				//读入一批消元子
				for (int j = 0; j < eliminatorWindowSize; j++)
				{
					//读取每一行稀疏向量形式的消元子
					string eliminatorSparseLine;
					if (eliminator.eof())
					{
						eliminator.close(); //这里使用seekg的方法调整读入位置会读取空行，因此重新打开
						eliminator.open(eliminatorPath, ios::in);
					}
					getline(eliminator, eliminatorSparseLine);

					//将读取到的压入消元子滑动窗中
					if (eliminatorSparseLine != "\0")
					{
						eliminatorSparseWindow.push_back(eliminatorSparseLine);
					}
				}

				//将字符串型消元子滑动窗转化为BitMap形式
				for (int j = 0; j < eliminatorSparseWindow.size(); j++)
				{
					BitMap eliminatorBitMap;
					eliminatorBitMap.init();
					eliminatorBitMap.sparseRowToBitMap(eliminatorSparseWindow[j]);
					eliminatorWindow.push_back(eliminatorBitMap);
					eliminatorWindowRange.insert(eliminatorBitMap.getHighestNumber());
				}

				int originSize = eliminatorWindow.size();

				//清空字符串型消元子滑动窗，节省内存
				vector<string>().swap(eliminatorSparseWindow);

				/*cout << "current "<<eliminatantWindow[i].toString() << endl;*/

				/*cout << "set ";
				unordered_set<int>::iterator it = eliminatorLeftHashSet.begin();
				while (it != eliminatorLeftHashSet.end())
				{
					cout << *it << ' ';
					it++;
				}
				cout << endl;*/

				/*cout << "range ";
				unordered_set<int>::iterator it1 = windowRange.begin();
				while (it1 != windowRange.end())
				{
					cout << *it1 << ' ';
					it1++;
				}
				cout << endl;*/

				//先判断当前行是否可以进行消元
				while (eliminatorWindowRange.find(eliminatantWindow[i].getHighestNumber())!=eliminatorWindowRange.end())
				{
					for (int j = 0; j < eliminatorWindow.size(); j++)
					{
						if (eliminatorWindow[j].getHighestNumber() == eliminatantWindow[i].getHighestNumber())
						{
							// eliminatantWindow[i].xorNormal(eliminatorWindow[j]);
							eliminatantWindow[i].xorSIMD(eliminatorWindow[j]);
							break;
						}
					}
				}

				/*string t1;
				eliminatantWindow[i].toString(t1);
				cout << "after:" << t1 << endl;*/

				//cout << "originSize" << originSize << endl;
				
				//判断当前行是否可以升格
				if (eliminatantWindow[i].getHighestNumber()!=-1 && eliminatorLeftHashSet.find(eliminatantWindow[i].getHighestNumber()) == eliminatorLeftHashSet.end())
				{
					//可以升格，加入到消元子中，同时更新消元子索引，同时表示该行运算完毕，写入结果文件
					eliminatantWindow[i].setToEliminator();
					eliminatorWindow.push_back(eliminatantWindow[i]);
					eliminatorLeftHashSet.insert(eliminatantWindow[i].getHighestNumber());
					writeEliminatorWindowToTempFile(eliminatorWindow,originSize);
					writeEliminatantWindowToResultFile(eliminatantWindow[i]);
				}

				//如果是空行，写入到结果文件中
				if (eliminatantWindow[i].getHighestNumber() == -1)
				{
					writeEliminatantWindowToResultFile(eliminatantWindow[i]);
				}
				
				vector<BitMap>().swap(eliminatorWindow);
				eliminatorWindowRange.clear();

				/*string temp;
				eliminatantWindow[i].toString(temp);
				cout << temp << endl;*/
			}
		}
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
	void writeEliminatorWindowToTempFile(vector<BitMap>& eliminatorWindow, int originSize)
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

	int i1 = exampleDirectory.find("矩阵列数") + 8;
	int i2 = exampleDirectory.find("，非零消元子");
	int len1 = i2 - i1;
	string columns = exampleDirectory.substr(i1,len1);
	stringstream ss;
	ss << columns;
	ss >> maxColunmAmount;

	GrobnerBasedGaussElimination g;
	g.init(exampleDirectoryPath);

	return 0;
}
