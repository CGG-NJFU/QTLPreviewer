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

#include "IMQTLPreviewer.h"

using namespace std;

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
 * 打印样本的表现型数据
 * @param data 样本表现型数组
 * @param size 样本大小
 * @param u0 统计量mu
 * @param s0 统计量sigma平方
 */
void printExpressData(EXPData *data, int size, EXPData u0, EXPData s0) {
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
 * @param u0 统计量mu 即平均值
 * @param s0 统计量sigma平方 及方差
 * @param ifPrintLog 是否打印日志
 */
void initExpressData(const string fileName, const int sampleNumber,
		EXPData* data, EXPData* u0, EXPData* s0, int ifPrintLog = VERBOSE_MODE) {
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
 * @param data 样本基因型数据
 * @param sampleNumber 样本大小
 * @param traitNumber 位点数量
 */
void printGeneData(string* data, int sampleNumber, int traitNumber) {
	cout << "=======MarkData======" << endl;
	int i, j;
	for (i = 0; i < traitNumber; i++) {
		cout << "[" << i << "]\t";
		for (j = 0; j < sampleNumber; j++) {
			if (j % 20 == 0)
				cout << endl << "[" << j << "]";
			cout << data[i][j];
		}
		cout << endl;
	}
}

/**
 * 初始化样本的基因型数据
 * @param fileName 样本数据文件的文件名
 * @param sampleNumber 样本大小
 * @param traitNumber 位点数量
 * @param data 数据数组
 * @param ifPrintLog 是否打印日志
 */
void initGeneData(string fileName, const int sampleNumber,
		const int traitNumber, string* data, int ifPrintLog = VERBOSE_MODE) {
	fstream fin;
	fin.open(fileName.data(), ios::in);

	for (int i = 0; i < traitNumber; i++) {
		string line;
		getline(fin, line);

		for (int j = 0; j < sampleNumber; j++) {
			data[i][j] = line[j];
		}
	}
	fin.close();

	if (ifPrintLog) {
		printGeneData(data, sampleNumber, traitNumber);
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
 * 将HAB型的基因数据转成ab型的基因数据，H-ab，A-aa，B-bb
 * @param input 输入数据
 * @return 转换后的数据
 */
string HAB2ab(const char input) {
	switch(input) {
		case 'H': return "ab";
		case 'A': return "aa";
		case 'B': return "bb";
		default: cout<<input; return NULL;
	}
}

/**
 * 查找基因型的条件概率
 * @param geneMatrix 基因型
 * @param rMatrix 基因型概率
 * @param find 查找的匹配基因型
 * @param ifPrintLog 是否打印日志
 * @return 基因型的条件概率
 */
double findGeneCP(const string* geneMatrix, const double* rMatrix, const string find, int ifPrintLog = VERBOSE_MODE) {
	double same=0; //该基因型概率之和
	double all=0; //总概率之和，相当书上表格横行之和
	for (int i=0;i<GENE_CP_SIZE;i++) {
		if ( geneMatrix[i].substr(0,GENE_CP_CUT) == find.substr(0,GENE_CP_CUT) ) {
			if (ifPrintLog) {
				cout <<geneMatrix[i] <<" added:\t" <<rMatrix[i] <<endl;
			}
			all += rMatrix[i];

			if (geneMatrix[i]==find) {
				same += rMatrix[i];
			}
		}
	}

	double re;
	if ( all == 0 ) re = 0; else re = same/all;

	if (ifPrintLog) {
		cout <<"CP: " <<re <<" (" <<same <<"/" <<all <<")" <<endl;
	}
	return re;
}

/**
 * 根据基因型数据统计对应的实际分布概率
 * @param _gp 实际分布概率
 * @param mk 基因型数据数组
 * @param sampleSize 样本大小
 * @param sampleIndex1 样本编号1
 * @param sampleIndex2 样本编号2
 * @param geneMatrix 基因型
 * @param rMatrix 条件型概率
 * @param ifUseShortGeneData 是否使用短基因表示
 * @param ifPrintLog 是否打印日志
 */
void calculateGP(double** _gp, string* mk, const int sampleSize, const int sampleIndex1,
		const int sampleIndex2, const string* geneMatrix, const double* rMatrix, bool ifUseShortGeneData, int ifPrintLog = VERBOSE_MODE) {
	double* gp = (double *) _gp;

	for (int i = 0; i < sampleSize; i++) {
		//亲本的基因型
		//TODO 修正，这里应该是string
		char cSample1 = mk[sampleIndex1][i];
		char cSample2 = mk[sampleIndex2][i];

		string gene;
		if (ifUseShortGeneData) {
			gene = HAB2ab(cSample1)+HAB2ab(cSample2);
		} else {
			gene = string(cSample1,1)+string(cSample2,1);
		}
		string qtl[3] = {"QQ", "Qq", "qq" };

		for (int j=0; j<3; j++) {
			*(gp + i * sampleSize + j) = findGeneCP(geneMatrix, rMatrix, gene+qtl[j], ifPrintLog);
		}
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
 * @param expData 表现型数据
 * @param s0 表现型的统计量sigima平方
 * @param s1
 * @param u0 表现型的统计量mu
 * @param u1
 * @param u2
 * @param u3
 * @param ifPrintLog 是否打印日志
 * @return LOD计分
 */
double calculateLOD(double** _gp, int sampleSize, EXPData* expData, EXPData s0,
		EXPData s1, EXPData u0, EXPData u1, EXPData u2, EXPData u3, int ifPrintLog =
				VERBOSE_MODE) {
	double* gp = (double*) _gp;
	double sum1 = 0, sum2 = 0;
	//cout <<"u0=" <<u0 <<endl <<"s0=" <<s0 <<endl;
	for (int i = 0; i < sampleSize; i++) {
		sum1 += log10(f(expData[i], u0, s0));
		sum2 += log10(
				*(gp + i * sampleSize + 0) * f(expData[i], u1, s1)
						+ *(gp + i * sampleSize + 1) * f(expData[i], u2, s1)
						+ *(gp + i * sampleSize + 2) * f(expData[i], u3, s1));
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
 * @param expData 表现型数据
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
void EMCalculate(EXPData* expData, double** _gp, const int sampleSize, EXPData u0,
		EXPData s0, EXPData* u1, EXPData* u2, EXPData* u3, EXPData* s1,
		int ifPrintLog = VERBOSE_MODE) {
	double *gp = (double*) _gp;

	// u1, u2, u3, sigma2 见113
	EXPData u10 = u0;
	EXPData u20 = u0;
	EXPData u30 = u0;
	EXPData s10 = s0;	//*s1+1;

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
			double gpi1f = *(gp + i * sampleSize + 0) * f(expData[i], u10, s10);
			double gpi2f = *(gp + i * sampleSize + 1) * f(expData[i], u20, s10);
			double gpi3f = *(gp + i * sampleSize + 2) * f(expData[i], u30, s10);
			double gpi = gpi1f + gpi2f + gpi3f;

			py1 += gpi1f / gpi;	//PI1的值
			py2 += gpi2f / gpi;	//PI2的值
			py3 += gpi3f / gpi;	//PI3的值

			py10 += gpi1f / gpi * expData[i];	//PI1和PI1*YI的和
			py20 += gpi2f / gpi * expData[i];	//PI2和PI2*YI的和
			py30 += gpi3f / gpi * expData[i];	//PI3和PI3*YI的和

			py4 += gpi1f / gpi * (expData[i] - u10) * (expData[i] - u10);
			py5 += gpi2f / gpi * (expData[i] - u20) * (expData[i] - u20);
			py6 += gpi3f / gpi * (expData[i] - u30) * (expData[i] - u30);
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
 * 计算基因型概率
 * @param gene 基因型的二进制表示，其中1表示大写，0表示小写
 * @param r1 左侧重组率
 * @param r2 右侧重组率
 * @param r 整体重组率，默认为0，通过r1和r2计算
 * @param zeta 符合系数，默认为1
 * @return 基因型概率
 */
double calcCross(bitset<GENE_CP_BIT> geneBit, double r1, double r2, double r=0, double zeta=1) {
	if (geneBit.size()!=GENE_CP_BIT) return 0;

	double r12 = r1*r2*zeta;
	if (r==0) {
		r = r1+r2-2*r12;
	}

	// 配子产生的概率，见书105页
	switch(geneBit.to_ulong()) {
		case 0: //M1QM2 000
		case 7: //m1qm2 111
			return (1.0-r1-r2+r12);// `=(1-r1)(1-r2);
		case 1: //M1Qm2 001
		case 6: //m1qM2 110
			return (r2-r12);// `=r2(1-r1);
		case 4: //M1qm2 100
		case 3: //m1QM2 011
			return (r1-r12);// `=r1(1-r2);
		case 5: //M1qM2 101
		case 2: //m1Qm2 010
			return r12;// `=r1*r2;
	}

	return 0;
}

/**
 * 打印基因型条件概率分布全数据
 * @param geneMatrix 基因型
 * @param rMatrix 条件概率
 */
void printGeneCP(string geneMatrix[GENE_CP_SIZE], double rMatrix[GENE_CP_SIZE]) {
	cout <<"1122:12 ->\t1122:12\tConditional Probability\n";
	for (int i=0; i<GENE_CP_SIZE; i++) {
		bitset<4> bs_up(i/GENE_CP_CUT);
		bitset<2> bs_down(i%GENE_CP_CUT);
		cout <<bs_up <<":" <<bs_down <<" -> \t"
				<<geneMatrix[i].substr(0,GENE_CP_CUT) <<":" << geneMatrix[i].substr(GENE_CP_CUT,2) <<"\t"
				<<rMatrix[i] <<endl;
	}
}

/**
 * 重排一个字符串内字符，如bdcae重排为abcde
 * @param input 需要重排的字符串
 * @return
 */
string reorderStr(string input) {
	sort(input.begin(), input.end());
	return input;
}

/**
 * 重排基因位顺序，如ba被重排为ab
 * @param geneMatrix 基因型数组
 * @return 成功运行返回1
 */
int reorderGene(string geneMatrix[GENE_CP_SIZE]) {
	for (int i=0; i<GENE_CP_SIZE; i++) {
			geneMatrix[i] = reorderStr(geneMatrix[i].substr(0,2)) +
					reorderStr(geneMatrix[i].substr(2,2)) +
					reorderStr(geneMatrix[i].substr(4,2));
	}
	return 1;
}

/**
 * 计算基因型条件概率分布全数据
 * @param geneMatrix 基因型
 * @param rMatrix 条件型概率
 * @param f 父本基因型
 * @param m 母本基因型
 * @param q QTL基因型
 * @param r1 左侧重组率
 * @param r2 右侧重组率
 * @param r 整体重组率
 * @param zeta 符合系数
 * @param ifReorder 是否自动重排基因位顺序，如ba被重排为ab
 * @param ifPrintLog 是否打印日志
 * @return 成功运行则返回1
 */
int calcGeneCP(string geneMatrix[GENE_CP_SIZE], double rMatrix[GENE_CP_SIZE], const string f, const string m, const string q, double r1, double r2, double r=0, double zeta=1, bool ifReorder=true, int ifPrintLog = VERBOSE_MODE) {
	for (int i=0; i<GENE_CP_SIZE; i++) {
		bitset<4> bs_fm(i/GENE_CP_CUT); //前4位为亲本基因型
		bitset<2> bs_qtl(i%GENE_CP_CUT); //后2位为QTL基因型

		int fGene = bs_fm[3]*4 + bs_qtl[1]*2 + bs_fm[1]; //父本基因型
		int qGene = bs_fm[2]*4 + bs_qtl[0]*2 + bs_fm[0]; //母本基因型

		bitset<3> bs_up(fGene);
		bitset<3> bs_down(qGene);

		rMatrix[i] = calcCross(bs_up, r1, r2, r, zeta) * calcCross(bs_down, r1, r2, r, zeta);

		geneMatrix[i] =
				(i&0x20?f.substr(1,1):f.substr(0,1)) +
				(i&0x10?f.substr(3,1):f.substr(2,1)) +
				(i&0x08?m.substr(1,1):m.substr(0,1)) +
				(i&0x04?m.substr(3,1):m.substr(2,1)) +
				(i&0x02?q.substr(1,1):q.substr(0,1)) +
				(i&0x01?q.substr(1,1):q.substr(0,1));
	}

	if (ifReorder) {
		reorderGene(geneMatrix);
	}

	if (ifPrintLog) {
		printGeneCP(geneMatrix, rMatrix);
	}
	return 1;
}

/**
 * 区间QTL分析
 * @param currentTrait 当前位点序号
 * @param u0 统计量mu
 * @param s0 统计量sigma平方
 * @param length 位点间隔
 * @param startPoint 开始位置
 * @param _mk 位点基因型数据
 * @param expData 位点表现型数据
 * @param sampleSize 样本大小
 * @param f 父本基因型
 * @param m 母本基因型
 * @param q QTL基因型
 * @param ifUseShortGeneData 是否使用短基因表示
 * @param ifPrintLog 是否打印日志
 * @return 返回LOD值
 */
double intervalQTL(int currentTrait, EXPData u0, EXPData s0, double length,
		double startPoint, string* mk, EXPData* expData, int sampleSize,
		string f, string m, string q, bool ifUseShortGeneData=false,
		int ifPrintLog = VERBOSE_MODE) {
	double r = d2r(length);
	double r1 = d2r(startPoint);
	double r2 = d2r(length - startPoint);

	if (ifPrintLog) {
		cout << "r=" << r << "\tr1=" << r1 << "\tr2=" << r2 << "\td1="
				<< startPoint << "\td2=" << (length - startPoint) << "\td="
				<< length << endl;
	}

	//不再使用该方法
	//无干扰时F2群体 QTL 基因频率分布数组函数p 书107页//
	//double p[9][3];
	//calculatePij(p, r, r1, r2, ifPrintLog);

	//计算基因型条件概率分布全数据，相当于生成书107页表
	double* rMatrix = new double[GENE_CP_SIZE];
	string* geneMatrix = new string[GENE_CP_SIZE];
	calcGeneCP(geneMatrix, rMatrix, f, m, q, r1, r2, r, 1, true, ifPrintLog);

	//统计各类型分布比率
	double **gp = new double*[sampleSize];
	for (int i = 0; i < sampleSize; i++) {
		gp[i] = new double[3];
	}

	calculateGP(gp, mk, sampleSize, currentTrait, currentTrait + 1, geneMatrix, rMatrix, ifUseShortGeneData, ifPrintLog);

	double u1 = 0, u2 = 0, u3 = 0, s1 = 0;
	// EM算法
	EMCalculate(expData, gp, sampleSize, u0, s0, &u1, &u2, &u3, &s1, ifPrintLog);

	return calculateLOD(gp, sampleSize, expData, s0, s1, u0, u1, u2, u3, ifPrintLog);
}

/**
 * 打印帮助信息
 * @param fileName 文件名
 */
void printHelp(char* fileName) {
	cout << "Usage:\t" << fileName << " [ConfigFileName]" << endl;
}

int mainQTL(int args, char* argv[]) {
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
		inputValueFromKeyboard("the Number of Samples", &iSampleSize);
		inputValueFromKeyboard("the Number of Traits", &iTraitNumber);
		inputValueFromKeyboard("the length of step", &step);
		inputValueFromKeyboard("the Filename of Express Data", &sExpressDataFile);
		inputValueFromKeyboard("the Filename of Gene Data", &sGeneDataFile);
		inputValueFromKeyboard("the Filename of Trait Interval Data", &sTraitIntervalFile);

		ifPrintInitDataReport = inputBooleanFromKeyboard("Do you want detailed Initialization Report?");
		ifPrintCalReport = inputBooleanFromKeyboard("Do you want detailed Calculation Report?");
		ifPrintFinalReport = inputBooleanFromKeyboard("Do you want detailed Final Report?");
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
			<< " traits in following files:" << endl
			<< "Express Data in:        " << sExpressDataFile << endl
			<< "Gene Data in:           " << sGeneDataFile << endl
			<< "Trait Interval Data in: " << sTraitIntervalFile << endl
			<< "Log Detail: " << endl << "\tInitialization Report - "
			<< (ifPrintInitDataReport ? "Yes" : "No") << endl
			<< "\tCalculation Report - " << (ifPrintCalReport ? "Yes" : "No")
			<< endl << "\tFinal Report - "
			<< (ifPrintFinalReport ? "Yes" : "No") << endl
			<< "Calculate step:\t" << step << endl;

	/// 读取表型数据文件，并计算表型数组的均值、方差
	EXPData* expData = new EXPData[iSampleSize];
	EXPData u0 = 0, s0 = 0;
	initExpressData(sExpressDataFile, iSampleSize, expData, &u0, &s0,
			ifPrintInitDataReport);

	/// 初始化基因型数据
	string* mk = new string[iTraitNumber];
	for (int i = 0; i < iTraitNumber; i++) {
		mk[i] = string(iSampleSize,' ');
	}
	/// 读取基因型数据文件
	initGeneData(sGeneDataFile, iSampleSize, iTraitNumber, mk, ifPrintInitDataReport);

	/// 初始化并读取位点距离数据，来源为见书115页
	double* traitInterval = new double[iTraitNumber];
	initIntervalData(sTraitIntervalFile, iTraitNumber - 1, traitInterval,
			ifPrintInitDataReport);

	/// 初始化亲本的基因型
	string fGene = "abab";
	string mGene = "abab";
	string qGene = "Qq";
	bool ifUseShortGeneData = true;

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
					traitInterval[currentTrait], startPoint, mk, expData, iSampleSize, fGene, mGene, qGene, ifUseShortGeneData ,ifPrintCalReport);
			if (ifPrintFinalReport) {
				cout << "[" << startPoint << "~" << (startPoint + step) << "]:["
						<< traitInterval[currentTrait] << "] LOD=" << LOD
						<< endl;
			}
		}
	}
	return 1;
}

int main(int args, char* argv[]) {
	getTimeStamp();
	int re = mainQTL(args, argv);
	float t = getTimeStamp();

	cout <<"Total time: " <<t <<" ms" <<endl
		<<"Average time:" << t <<" ms" <<endl;

	return re;
}
