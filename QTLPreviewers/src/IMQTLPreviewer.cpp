//============================================================================
// Name        : IMQTLPreviewer.cpp
// Author      : Ethan
// Version     : Master
// Copyright   : Ethan @2012
//		iSampleSize = 189;
//		iTraitNumber = 9;
//		step=1.0;
//
//		sExpressDataFile = "data/express-66.txt";
//		sGeneDataFile = "data/gene-66.txt";
//		sTraitIntervalFile = "data/traitInterval-66.txt";
//
//		ifPrintInitDataReport = false;
//		ifPrintCalReport = false;
//		ifPrintFinalReport = true;
// Description : Previewer for Interval Mapping of QTL Analysis
//============================================================================

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define ACCURACY_REQ 0.000001
#define VERBOSE_MODE true

/**
 * 概率密度函数
 * @param x 输入
 * @param u mu
 * @param s2 sigma平方
 * @return
 */
double f(double x, double u, double s2) {
	return ((1 / sqrt(2 * M_PI * s2)) * exp(-0.5 * (x - u) * (x - u) / s2));
}

/**
 * 获取平均值
 * @param data 浮点数组
 * @param size 数组尺寸
 * @return 数组数据的平均值
 */
double getAverage(double* data, int size) {
	int i;
	double sum = 0;
	for (i = 0; i < size; i++)
		sum += data[i];
	return sum / size;
}

/**
 * 获取方差
 * @param data 浮点数组
 * @param size 数组尺寸
 * @param average 平均值，默认为0，视为自动计算
 * @return 数组数据的方差
 */
double getVariance(double* data, int size, double average = 0) {
	int i;
	double sum = 0;
	if (0 == average) {
		average = getAverage(data, size);
	}
	for (i = 0; i < size; i++) {
		sum += (data[i] - average) * (data[i] - average);
	}
	return sum / size;
}

/**
 * 打印样本的表现型数据
 * @param data 样本表现型数组
 * @param size 样本大小
 * @param u0 统计量mu
 * @param s0 统计量sigma平方
 */
void printExpressData(double *data, int size, double u0, double s0) {
	cout << "=======ExpressData=======" << endl;
	for (int i = 0; i < size; i++) {
		cout << "[" << i << "]\t" << data[i] << endl;
	}
	cout << "u0=" << u0 << endl << "s0=" << s0 << endl;
}

/**
 * 初始化样本的表现型数据
 * @param fileName 样本数据文件的文件名
 * @param sampleNumber 样本大小
 * @param data 数据数组
 * @param u0 统计量mu
 * @param s0 统计量sigma平方
 * @param ifPrintLog 是否打印日志
 */
void initExpressData(const string fileName, const int sampleNumber,
		double* data, double* u0, double* s0, int ifPrintLog = VERBOSE_MODE) {
	fstream fin;
	fin.open(fileName.data(), ios::in);
	for (int i = 0; i < sampleNumber; i++) {
		fin >> data[i];
	}
	fin.close();

	*u0 = getAverage(data, sampleNumber);
	*s0 = getVariance(data, sampleNumber, *u0);
	if (ifPrintLog) {
		printExpressData(data, sampleNumber, *u0, *s0);
	}
}

/**
 * 打印样本的基因型数据
 * @param _data 样本基因型数据
 * @param sampleNumber 样本大小
 * @param traitNumber 位点数量
 */
void printGeneData(char** _data, int sampleNumber, int traitNumber) {
	char* data = (char*) _data;
	cout << "=======MarkData======" << endl;
	int i, j;
	for (i = 0; i < traitNumber; i++) {
		cout << "[" << i << "]\t";
		for (j = 0; j < sampleNumber; j++) {
			if (j % 20 == 0)
				cout << endl << "[" << j << "]";
			cout << *(data + i * sampleNumber + j);
		}
		cout << endl;
	}
}

/**
 * 初始化样本的基因型数据
 * @param fileName 样本数据文件的文件名
 * @param sampleNumber 样本大小
 * @param traitNumber 位点数量
 * @param _data 数据数组
 * @param ifPrintLog 是否打印日志
 */
void initGeneData(string fileName, const int sampleNumber,
		const int traitNumber, char** _data, int ifPrintLog = VERBOSE_MODE) {
	char* data = (char*) _data;
	fstream fin;
	fin.open(fileName.data(), ios::in);

	for (int i = 0; i < traitNumber; i++) {
		string line;
		getline(fin, line);

		for (int j = 0; j < sampleNumber; j++) {
			*(data + i * sampleNumber + j) = line[j];
		}
	}
	fin.close();

	if (ifPrintLog) {
		printGeneData(_data, sampleNumber, traitNumber);
	}
}

/**
 * 打印遗传图谱的遗传区间距离数据
 * @param data 区间距离数据
 * @param intervalNumber 区间个数
 */
void printIntervalData(double* data, int intervalNumber) {
	cout << "=======TraitIntervalData=======" << endl;
	for (int i = 0; i < intervalNumber; i++) {
		cout << "[" << i << "]\t" << data[i] << endl;
	}
}

/**
 * 初始化遗传图谱的遗传区间距离数据
 * @param fileName 区间距离数据文件的文件名
 * @param intervalNumber 区间个数
 * @param data 区间距离数据
 * @param ifPrintLog 是否打印日志
 */
void initIntervalData(string fileName, const int intervalNumber, double* data,
		int ifPrintLog = VERBOSE_MODE) {
	fstream fin;
	fin.open(fileName.data(), ios::in);
	for (int i = 0; i < intervalNumber; i++) {
		fin >> data[i];
	}
	fin.close();

	if (ifPrintLog) {
		printIntervalData(data, intervalNumber);
	}
}

/**
 * 无干扰时F2群体QTL基因频率分布数组函数p，见书107页
 * @param p 基因频率分布数组函数
 * @param r 整体遗传距离
 * @param r1 遗传距离1
 * @param r2 遗传距离2
 * @param ifPrintLog 是否打印日志
 */
void calculatePij(double p[9][3], double r, double r1, double r2,
		int ifPrintLog = VERBOSE_MODE) {
	p[0][0] = (1 - r1) * (1 - r1) * (1 - r2) * (1 - r2) / ((1 - r) * (1 - r));
	p[0][1] = 2 * r1 * r2 * (1 - r1) * (1 - r2) / ((1 - r) * (1 - r));
	p[0][2] = r1 * r1 * r2 * r2 / ((1 - r) * (1 - r));
	p[1][0] = (1 - r1) * (1 - r1) * r2 * (1 - r2) / (r * (1 - r));
	p[1][1] = r1 * (1 - r1) * (1 - 2 * r2 + 2 * r2 * r2) / (r * (1 - r));
	p[1][2] = r1 * r1 * r2 * (1 - r2) / (r * (1 - r));
	p[2][0] = (1 - r1) * (1 - r1) * r2 * r2 / (r * r);
	p[2][1] = 2 * r1 * r2 * (1 - r1) * (1 - r2) / (r * r);
	p[2][2] = r1 * r1 * (1 - r2) * (1 - r2) / (r * r);
	p[3][0] = r1 * (1 - r1) * (1 - r2) * (1 - r2) / (r * (1 - r));
	p[3][1] = (1 - 2 * r1 + 2 * r1 * r1) * r2 * (1 - r2) / (r * (1 - r));
	p[3][2] = r1 * (1 - r1) * r2 * r2 / (r * (1 - r));
	p[4][0] = 2 * r1 * r2 * (1 - r1) * (1 - r2) / (1 - 2 * r + 2 * r * r);
	p[4][1] = (1 - 2 * r1 + 2 * r1 * r1) * (1 - 2 * r2 + 2 * r2 * r2)
			/ (1 - 2 * r + 2 * r * r);
	p[4][2] = p[4][0];
	p[5][0] = p[3][2];
	p[5][1] = p[3][1];
	p[5][2] = p[3][0];
	p[6][0] = p[2][2];
	p[6][1] = p[2][1];
	p[6][2] = p[2][0];
	p[7][0] = p[1][2];
	p[7][1] = p[1][1];
	p[7][2] = p[1][0];
	p[8][0] = p[0][2];
	p[8][1] = p[0][1];
	p[8][2] = p[0][0];

	if (ifPrintLog) {
		for (int i = 0; i < 9; i++) {
			double sum = 0;
			for (int j = 0; j < 3; j++) {
				sum += p[i][j];
				cout << p[i][j] << "\t";
			}
			cout << sum << "=1.0" << endl;
		}
	}
}

/**
 * 根据基因型数据统计对应的实际分布概率
 * @param _gp 实际分布概率
 * @param _mk 基因型数据数组
 * @param sampleSize 样本大小
 * @param sampleIndex1 样本编号1
 * @param sampleIndex2 样本编号2
 * @param p 分布概率函数
 * @param ifPrintLog 是否打印日志
 */
void calculateGP(double** _gp, char** _mk, int sampleSize, int sampleIndex1,
		int sampleIndex2, double p[9][3], int ifPrintLog = VERBOSE_MODE) {
	double* gp = (double *) _gp;
	char* mk = (char*) _mk;

	for (int i = 0; i < sampleSize; i++) {
		int lineNumber = 0;

		char cSample1 = *(mk + sampleIndex1 * sampleSize + i);
		char cSample2 = *(mk + sampleIndex2 * sampleSize + i);

		if (cSample1 == 'A') {
			if (cSample2 == 'A')
				lineNumber = 0;
			else if (cSample2 == 'H')
				lineNumber = 1;
			else if (cSample2 == 'B')
				lineNumber = 2;
		} else if (cSample1 == 'H') {
			if (cSample2 == 'A')
				lineNumber = 3;
			else if (cSample2 == 'H')
				lineNumber = 4;
			else if (cSample2 == 'B')
				lineNumber = 5;
		} else if (cSample1 == 'B') {
			if (cSample2 == 'A')
				lineNumber = 6;
			else if (cSample2 == 'H')
				lineNumber = 7;
			else if (cSample2 == 'B')
				lineNumber = 8;
		}

		*(gp + i * sampleSize + 0) = p[lineNumber][0];
		*(gp + i * sampleSize + 1) = p[lineNumber][1];
		*(gp + i * sampleSize + 2) = p[lineNumber][2];
		//cout<<i <<"\t" <<cSample1 <<"\t" <<cSample2 <<"\t" <<lineNumber <<endl;
	}

	if (ifPrintLog) {
		cout << "=======GPData======" << endl;
		for (int c = 0; c < sampleSize; c++) {
			cout << c << "\t";
			for (int l = 0; l < 3; l++)
				cout << *(gp + c * sampleSize + l) << "\t";
			cout << endl;
		}
	}
}

/**
 * 利用LOD计分法计算LOD值
 * @param _gp 实际分布概率
 * @param sampleSize 样本大小
 * @param y 表现型数据
 * @param s0 表现型的统计量sigima平方
 * @param s1
 * @param u0 表现型的统计量mu
 * @param u1
 * @param u2
 * @param u3
 * @param ifPrintLog 是否打印日志
 * @return LOD计分
 */
double calculateLOD(double** _gp, int sampleSize, double* y, double s0,
		double s1, double u0, double u1, double u2, double u3, int ifPrintLog =
				VERBOSE_MODE) {
	double* gp = (double*) _gp;
	double sum1 = 0, sum2 = 0;
	//cout <<"u0=" <<u0 <<endl <<"s0=" <<s0 <<endl;
	for (int i = 0; i < sampleSize; i++) {
		sum1 += log10(f(y[i], u0, s0));
		sum2 += log10(
				*(gp + i * sampleSize + 0) * f(y[i], u1, s1)
						+ *(gp + i * sampleSize + 1) * f(y[i], u2, s1)
						+ *(gp + i * sampleSize + 2) * f(y[i], u3, s1));
	}

	// LOD=log(L0/L)=logL0-LogL
	// 连乘Log改为加法，除法改减法
	double lod;
	lod = sum2 - sum1;

	if (ifPrintLog) {
		cout << "sum1=" << sum1 << "\tsum2=" << sum2 << "\tLOD=" << lod << endl;
	}
	return lod;
}

/**
 * 利用EM算法迭代计算最大似然估计
 * @param y 表现型数据
 * @param _gp 实际分布概率
 * @param sampleSize 样本大小
 * @param u0 样本的统计量mu
 * @param s0 样本的统计量sigma平方
 * @param u1
 * @param u2
 * @param u3
 * @param s1
 * @param ifPrintLog 是否打印日志
 */
void EMCalculate(double* y, double** _gp, const int sampleSize, double u0,
		double s0, double* u1, double* u2, double* u3, double* s1,
		int ifPrintLog = VERBOSE_MODE) {
	double *gp = (double*) _gp;

	// u1, u2, u3, sigma2 见113
	double u10 = u0;
	double u20 = u0;
	double u30 = u0;
	double s10 = s0;	//*s1+1;

	// 收敛精度要求
	int step = 1;
	if (ifPrintLog) {
		cout << "=======EMSteps======" << endl;
		cout << "[step]\tu1\tu2\tu3\ts1\n";
	}
	while (fabs(u10 - *u1) > ACCURACY_REQ || fabs(u20 - *u2) > ACCURACY_REQ
			|| fabs(u30 - *u3) > ACCURACY_REQ || fabs(s10 - *s1) > ACCURACY_REQ) {
		double py1 = 0, py10 = 0;
		double py2 = 0, py20 = 0;
		double py3 = 0, py30 = 0;
		double py4 = 0, py5 = 0, py6 = 0;

		int i;
		for (i = 0; i < sampleSize; i++) {
			double gpi1f = *(gp + i * sampleSize + 0) * f(y[i], u10, s10);
			double gpi2f = *(gp + i * sampleSize + 1) * f(y[i], u20, s10);
			double gpi3f = *(gp + i * sampleSize + 2) * f(y[i], u30, s10);
			double gpi = gpi1f + gpi2f + gpi3f;

			py1 += gpi1f / gpi;		//PI1的值
			py2 += gpi2f / gpi;	//PI2的值
			py3 += gpi3f / gpi;	//PI3的值

			py10 += gpi1f / gpi * y[i];	//PI1和PI1*YI的和
			py20 += gpi2f / gpi * y[i];	//PI2和PI2*YI的和
			py30 += gpi3f / gpi * y[i];	//PI3和PI3*YI的和

			py4 += gpi1f / gpi * (y[i] - u10) * (y[i] - u10);
			py5 += gpi2f / gpi * (y[i] - u20) * (y[i] - u20);
			py6 += gpi3f / gpi * (y[i] - u30) * (y[i] - u30);
		}

		*u1 = u10;
		u10 = py10 / py1;
		*u2 = u20;
		u20 = py20 / py2;
		*u3 = u30;
		u30 = py30 / py3;
		*s1 = s10;
		s10 = (py4 + py5 + py6) / sampleSize;

		if (ifPrintLog) {
			cout << "[" << step++ << "]\t" << *u1 << "\t" << *u2 << "\t" << *u3
					<< "\t" << *s1 << endl;
		}
	}	//end-while

	//输出终值
	*u1 = u10;
	*u2 = u20;
	*u3 = u30;
	*s1 = s10;
	if (ifPrintLog) {
		cout << endl << "Finals:" << endl;
		cout << "u1=" << *u1 << "\tu2=" << *u2 << "\tu3=" << *u3 << "\ts1="
				<< *s1 << endl;
	}
}

/**
 * Haldane作图函数
 * @param d 遗传距离
 * @return 重组率
 */
double d2rHaldane(double d) {
	return 0.5 * (1 - exp(-2 * d * 0.01));
}

/**
 * Kosambi作图函数
 * @param d 遗传距离
 * @return 重组率
 */
double d2rKosambi(double d) {
	return 0.5 * tanh(2 * d * 0.01);
}

/**
 * 作图函数，默认使用Haldane作图函数
 * @param x 遗传距离
 * @return 重组率
 */
double d2r(double x) {
	return d2rHaldane(x);
}

/**
 * 区间QTL分析
 * @param currentTrait 当前位点序号
 * @param u0 统计量mu
 * @param s0 统计量sigma平方
 * @param length 位点间隔
 * @param startPoint 开始位置
 * @param _mk 位点基因型数据
 * @param y 位点表现型数据
 * @param sampleSize 样本大小
 * @param ifPrint 是否打印日志
 * @return
 */
double intervalQTL(int currentTrait, double u0, double s0, double length,
		double startPoint, char** _mk, double* y, int sampleSize,
		int ifPrintLog = VERBOSE_MODE) {
	double r = d2r(length);
	double r1 = d2r(startPoint);
	double r2 = d2r(length - startPoint);

	if (ifPrintLog) {
		cout << "r=" << r << "\tr1=" << r1 << "\tr2=" << r2 << "\td1="
				<< startPoint << "\td2=" << (length - startPoint) << "\td="
				<< length << endl;
	}

	//无干扰时F2群体 QTL 基因频率分布数组函数p 书107页//
	double p[9][3];
	calculatePij(p, r, r1, r2, ifPrintLog);

	//统计各类型分布比率
	double **gp = new double*[sampleSize];
	for (int i = 0; i < sampleSize; i++) {
		gp[i] = new double[3];
	}

	calculateGP(gp, _mk, sampleSize, currentTrait, currentTrait + 1, p,
			ifPrintLog);

	double u1 = 0, u2 = 0, u3 = 0, s1 = 0;
	// EM算法
	EMCalculate(y, gp, sampleSize, u0, s0, &u1, &u2, &u3, &s1, ifPrintLog);

	return calculateLOD(gp, sampleSize, y, s0, s1, u0, u1, u2, u3, ifPrintLog);
}

/**
 * 打印帮助信息
 * @param fileName 文件名
 */
void printHelp(char* fileName) {
	cout << "Usage:\t" << fileName << " [ConfigFileName]" << endl;
}

int main(int args, char* argv[]) {
	int iSampleSize;
	int iTraitNumber;
	double step;

	string sExpressDataFile;
	string sGeneDataFile;
	string sTraitIntervalFile;

	bool ifPrintInitDataReport;
	bool ifPrintCalReport;
	bool ifPrintFinalReport;

	if (args == 1) {
		/// 若参数个数为1，为手动模式运行，从操作台读取各个参数。
		cout << "Please input the Number of Samples:";
		cin >> iSampleSize;

		cout << iSampleSize << endl << "Please input the Number of Traits: ";
		cin >> iTraitNumber;

		cout << iTraitNumber << endl << "Please input the length of step: ";
		cin >> step;

		cout << step << endl << "Please locate the file of Express Data: ";
		cin >> sExpressDataFile;

		cout << sExpressDataFile << endl
				<< "Please locate the file of Gene Data: ";
		cin >> sGeneDataFile;

		cout << sGeneDataFile << endl
				<< "Please locate the file of Trait Interval Data: ";
		cin >> sTraitIntervalFile;

		cout << sTraitIntervalFile << endl
				<< "Do you want detailed Initialization Report?(0/1) ";
		cin >> ifPrintInitDataReport;

		cout << (ifPrintInitDataReport ? "True" : "False") << endl
				<< "Do you want detailed Calculation Report?(0/1) ";
		cin >> ifPrintCalReport;

		cout << (ifPrintCalReport ? "True" : "False") << endl
				<< "Do you want detailed Final Report?(0/1) ";
		cin >> ifPrintFinalReport;

		cout << (ifPrintFinalReport ? "True" : "False") << endl;
	} else if (args == 2) {
		/// 若参数个树为2，为自动模式，从配置文件读取参数。
		cout << "Input File:" << argv[1] << endl;
		fstream fin;
		fin.open(argv[1], ios::in);

		fin >> iSampleSize >> iTraitNumber >> step;
		fin >> sExpressDataFile;
		fin >> sGeneDataFile;
		fin >> sTraitIntervalFile;
		fin >> ifPrintInitDataReport >> ifPrintCalReport >> ifPrintFinalReport;

		fin.close();
	} else {
		/// 参数个数错误，打印帮助。
		printHelp(argv[0]);
	}

	/// 汇总显示各输入参数
	cout << "Input Summary:\t" << iSampleSize << " samples of " << iTraitNumber
			<< " traits in following files:" << endl;
	cout << "Express Data in:        " << sExpressDataFile << endl;
	cout << "Gene Data in:           " << sGeneDataFile << endl;
	cout << "Trait Interval Data in: " << sTraitIntervalFile << endl;
	cout << "Log Detail: " << endl << "\tInitialization Report - "
			<< (ifPrintInitDataReport ? "Yes" : "No") << endl
			<< "\tCalculation Report - " << (ifPrintCalReport ? "Yes" : "No")
			<< endl << "\tFinal Report - "
			<< (ifPrintFinalReport ? "Yes" : "No") << endl;
	cout << "Calculate step:\t" << step << endl;

	/// 读取表型数据文件，并计算表型数组的均值、方差
	double* y = new double[iSampleSize];
	double u0 = 0, s0 = 0;
	initExpressData(sExpressDataFile, iSampleSize, y, &u0, &s0,
			ifPrintInitDataReport);

	/// 初始化基因型数据
	char** mk = new char*[iTraitNumber];
	for (int i = 0; i < iTraitNumber; i++) {
		mk[i] = new char[iSampleSize];
	}
	/// 读取基因型数据文件
	initGeneData(sGeneDataFile, iSampleSize, iTraitNumber, mk,
			ifPrintInitDataReport);

	/// 初始化并读取位点距离数据，来源为见书115页
	double* traitInterval = new double[iTraitNumber];
	initIntervalData(sTraitIntervalFile, iTraitNumber - 1, traitInterval,
			ifPrintInitDataReport);

	/// 逐个位点计算LOD值
	for (int currentTrait = 0; currentTrait < iTraitNumber - 1;
			currentTrait++) {
		if (ifPrintFinalReport) {
			cout << "---------" << (currentTrait + 1) << "---("
					<< traitInterval[currentTrait] << ")--------" << endl;
		}
		/// 按步长推进计算LOD值
		for (double startPoint = 0.0; startPoint < traitInterval[currentTrait];
				startPoint += step) {
			double LOD = intervalQTL(currentTrait, u0, s0,
					traitInterval[currentTrait], startPoint, mk, y, iSampleSize, ifPrintCalReport);
			if (ifPrintFinalReport) {
				cout << "[" << startPoint << "~" << (startPoint + step) << "]:["
						<< traitInterval[currentTrait] << "] LOD=" << LOD
						<< endl;
			}
		}
	}
	return 1;
}
