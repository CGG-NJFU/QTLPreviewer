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

extern Category& logger;

using namespace std;
using namespace log4cpp;

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
 * 初始化样本的表现型数据
 * @param fileName 样本数据文件的文件名
 * @param sampleNumber 样本大小
 * @param expData 数据数组
 * @param u0 统计量mu 即平均值
 * @param s0 统计量sigma平方 及方差
 * @param ifPrintLog 是否打印日志
 */
void initExpressData(const string fileName, const int sampleNumber,
		vector<EXPData>& expData, EXPData* u0, EXPData* s0, const bool ifPrintLog = VERBOSE_MODE) {
	int realsize = readFile2Vector(fileName, expData);
	if ( sampleNumber != realsize ) {
		logger <<Priority::WARN <<"Sample size (" <<realsize <<") might be wrong, please check.";
	}

	*u0 = getAverage( expData.begin(), expData.end() );
	*s0 = getVariance( expData.begin(), expData.end(), *u0);

	if (ifPrintLog) {
		printExpressData(expData, *u0, *s0);
	}
}

/**
 * 初始化样本的基因型数据
 * @param fileName 样本数据文件的文件名
 * @param data 数据数组
 * @param ifPrintLog 是否打印日志
 */
void initChildrenGeneData(const string fileName, vector<vector<string> >& data, const bool ifPrintLog = VERBOSE_MODE) {
	readFile2Matrix(fileName, data);

	if (ifPrintLog) {
		printChildrenGeneData(data);
	}
}

/**
 * 初始化亲本的基因型数据
 * @param fileName 亲本数据文件的文件名
 * @param geneData 数据数组
 * @param traitNumber 位点数量
 * @param parentName 亲本名称
 * @param ifPrintLog 是否打印日志
 */
void initParentGeneData(const string fileName, vector<string>& geneData, const int traitNumber, const string parentName, const bool ifPrintLog = VERBOSE_MODE) {
	int realsize = readFile2Vector(fileName, geneData);
	if ( traitNumber != realsize ) {
		logger <<Priority::WARN <<"Trait number (" <<realsize <<") might be wrong, please check.";
	}

	if (ifPrintLog) {
		printParentGeneData(geneData, parentName);
	}
}

/**
 * 初始化遗传图谱的遗传区间距离数据
 * @param fileName 区间距离数据文件的文件名
 * @param intervalNumber 区间个数
 * @param data 区间距离数据
 * @param ifPrintLog 是否打印日志
 */
void initIntervalData(const string fileName, const int intervalNumber, vector<double>& data,
		const bool ifPrintLog = VERBOSE_MODE) {
	int realsize = readFile2Vector(fileName, data);
	if ( intervalNumber != realsize ) {
		logger <<Priority::WARN <<"The number of intervals (" <<realsize <<") might be wrong, please check.";
	}

	if (ifPrintLog) {
		printIntervalData(data);
	}
}

/**
 * 将HAB型的基因数据转成ab型的基因数据，H-ab，A-aa，B-bb
 * @param input 输入数据
 * @return 转换后的数据
 */
string HAB2ab(const string input) {
	if (input=="H") return "ab";
	else if (input=="A") return "aa";
	else if (input=="B") return "bb";
	else return NULL;
}

/**
 * 查找基因型的条件概率
 * @param geneMatrix 基因型
 * @param rMatrix 基因型概率
 * @param find 查找的匹配基因型
 * @param ifPrintLog 是否打印日志
 * @return 基因型的条件概率
 */
double findGeneCP(const vector<string> geneMatrix, const vector<double> rMatrix, const string find, const bool ifPrintLog = VERBOSE_MODE) {
	double same=0; //该基因型概率之和
	double all=0; //总概率之和，相当书上表格横行之和
	for (unsigned int i=0; i<geneMatrix.size(); i++) {
		if ( geneMatrix[i].substr(0,GENE_CP_CUT) == find.substr(0,GENE_CP_CUT) ) {
			if (ifPrintLog) {
				logger <<Priority::DEBUG <<geneMatrix[i] <<" added:\t" <<rMatrix[i];
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
		logger <<Priority::DEBUG <<"CP: " <<re <<" (" <<same <<"/" <<all <<")";
	}
	return re;
}

/**
 * 根据基因型数据统计对应的实际分布概率
 * @param gp 实际分布概率
 * @param mk 基因型数据数组
 * @param qtlGene QTL基因的字符串
 * @param sampleSize 样本大小
 * @param sampleIndex1 样本编号1
 * @param sampleIndex2 样本编号2
 * @param geneMatrix 基因型
 * @param rMatrix 条件型概率
 * @param ifPrintLog 是否打印日志
 */
void calculateGP(vector<vector<double> >& gp, const vector<vector<string> >& mk, const string qtlGene, const int sampleSize, const int sampleIndex1,
		const int sampleIndex2, const vector<string>& geneMatrix, const vector<double>& rMatrix, const bool ifPrintLog = VERBOSE_MODE) {
	for (int i = 0; i < sampleSize; i++) {
		vector<string> qtl(3);
		qtl[0] = qtlGene.substr(0,1)+qtlGene.substr(0,1); //"QQ"
		qtl[1] = qtlGene.substr(0,1)+qtlGene.substr(1,1); //"Qq"
		qtl[2] = qtlGene.substr(1,1)+qtlGene.substr(1,1); //"qq"

		///获取并拼接样本的基因型
		string gene = mk[sampleIndex1][i] + mk[sampleIndex2][i];

		///计算实际分布概率
		for (int j=0; j<3; j++) {
			gp[i][j] = findGeneCP(geneMatrix, rMatrix, gene+qtl[j], ifPrintLog);
		}
	}

	if (ifPrintLog) {
		logger <<Priority::DEBUG << "=======GPData======";
		for (int c = 0; c < sampleSize; c++) {
			logger <<Priority::DEBUG << c << "\t" <<gp[c][0] <<"\t" <<gp[c][1] <<"\t"<<gp[c][1];
		}
	}
}

/**
 * 利用LOD计分法计算LOD值
 * @param gp 实际分布概率
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
double calculateLOD(const vector<vector<double> > gp, const int sampleSize, const vector<EXPData> expData,
		const EXPData u0, const EXPData s0,
		const EXPData u1, const EXPData u2, const EXPData u3, const EXPData s1,
		const bool ifPrintLog = VERBOSE_MODE) {
	double sum1 = 0, sum2 = 0;
	logger <<Priority::DEBUG <<"u0=" <<u0 <<"\ts0=" <<s0;
	for (int i = 0; i < sampleSize; i++) {
		sum1 += log10(f(expData[i], u0, s0));
		sum2 += log10(
				gp[i][0] * f(expData[i], u1, s1)
						+ gp[i][1] * f(expData[i], u2, s1)
						+ gp[i][2] * f(expData[i], u3, s1));
	}

	// LOD=log(L0/L)=logL0-LogL
	// 连乘Log改为加法，除法改减法
	double lod;
	lod = sum2 - sum1;

	if (ifPrintLog) {
		logger <<Priority::DEBUG <<"sum1=" <<sum1 <<"\tsum2=" <<sum2 <<"\tLOD=" <<lod;
	}
	return lod;
}

/**
 * 利用EM算法迭代计算最大似然估计
 * @param expData 表现型数据
 * @param gp 实际分布概率
 * @param sampleSize 样本大小
 * @param u0 样本的统计量mu
 * @param s0 样本的统计量sigma平方
 * @param u1
 * @param u2
 * @param u3
 * @param s1
 * @param ifPrintLog 是否打印日志
 */
void EMCalculate(const vector<vector<double> >& gp, const int sampleSize, const vector<EXPData> expData,
		const EXPData u0, const EXPData s0, EXPData* u1, EXPData* u2, EXPData* u3, EXPData* s1,
		const bool ifPrintLog = VERBOSE_MODE) {
	// u1, u2, u3, sigma2 见113
	EXPData u10 = u0;
	EXPData u20 = u0;
	EXPData u30 = u0;
	EXPData s10 = s0;	//*s1+1;

	// 收敛精度要求
	int step = 1;
	if (ifPrintLog) {
		logger <<Priority::DEBUG <<"=======EMSteps======";
		logger <<Priority::DEBUG <<"[step]\tu1\tu2\tu3\ts1";
	}
	while (fabs(u10 - *u1) > ACCURACY_REQ || fabs(u20 - *u2) > ACCURACY_REQ
			|| fabs(u30 - *u3) > ACCURACY_REQ || fabs(s10 - *s1) > ACCURACY_REQ) {
		double py1 = 0, py10 = 0;
		double py2 = 0, py20 = 0;
		double py3 = 0, py30 = 0;
		double py4 = 0, py5 = 0, py6 = 0;

		int i;
		for (i = 0; i < sampleSize; i++) {
			double gpi1f = gp[i][0] * f(expData[i], u10, s10);
			double gpi2f = gp[i][1] * f(expData[i], u20, s10);
			double gpi3f = gp[i][2] * f(expData[i], u30, s10);
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
			logger <<Priority::DEBUG << "[" << step++ << "]\t" << *u1 << "\t" << *u2 << "\t" << *u3
					<< "\t" << *s1;
		}
	}	//end-while

	//输出终值
	*u1 = u10;
	*u2 = u20;
	*u3 = u30;
	*s1 = s10;
	if (ifPrintLog) {
		logger <<Priority::DEBUG << "Finals:\tu1=" <<*u1 <<"\tu2=" <<*u2 <<"\tu3=" <<*u3 <<"\ts1=" <<*s1;
	}
}

/**
 * 计算基因型概率
 * @param geneBit 基因型的二进制表示，其中1表示大写，0表示小写
 * @param r1 左侧重组率
 * @param r2 右侧重组率
 * @param r 整体重组率，默认为0，通过r1和r2计算
 * @param zeta 符合系数，默认为1
 * @return 基因型概率
 */
double calcCross(const bitset<GENE_CP_BIT> geneBit,
		const double r1, const double r2, const double r=0, const double zeta=1) {
	if (geneBit.size()!=GENE_CP_BIT) return 0;

	double r12 = r1*r2*zeta;
	double r0 = r;
	if (r==0) {
		r0 = r1+r2-2*r12;
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
 * 重排基因位顺序，如ba被重排为ab
 * @param geneMatrix 基因型数组
 * @return 成功运行返回1
 */
void reorderGene(vector<string>& geneMatrix) {
	for (unsigned int i=0; i<geneMatrix.size(); i++) {
			geneMatrix[i] = reorderStr(geneMatrix[i].substr(0,2)) +
					reorderStr(geneMatrix[i].substr(2,2)) +
					reorderStr(geneMatrix[i].substr(4,2));
	}
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
int calcGeneCP(vector<string>& geneMatrix, vector<double>& rMatrix,
		const string f, const string m, const string q,
		const double r1, const double r2, const double r=0, const double zeta=1,
		const bool ifReorder=true, const bool ifPrintLog = VERBOSE_MODE) {
	for (unsigned int i=0; i<geneMatrix.size(); i++) {
		bitset<4> bs_fm(i/GENE_CP_CUT); //前4位为亲本基因型
		bitset<2> bs_qtl(i%GENE_CP_CUT); //后2位为QTL基因型

		int fGene = bs_fm[3]*4 + bs_qtl[1]*2 + bs_fm[1]; //父本基因型
		int mGene = bs_fm[2]*4 + bs_qtl[0]*2 + bs_fm[0]; //母本基因型

		bitset<3> bs_up(fGene);
		bitset<3> bs_down(mGene);

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
 * @param mk 位点基因型数据
 * @param expData 位点表现型数据
 * @param sampleSize 样本大小
 * @param f 父本基因型
 * @param m 母本基因型
 * @param q QTL基因型
 * @param ifPrintLog 是否打印日志
 * @return 返回LOD值
 */
double intervalQTL(const int currentTrait, const EXPData u0, const EXPData s0, const double length,
		const double startPoint, const vector<vector<string> >& mk, const vector<EXPData> expData, const int sampleSize,
		const string f, const string m, const string q,
		const bool ifPrintLog = VERBOSE_MODE) {
	double r = d2r(length);
	double r1 = d2r(startPoint);
	double r2 = d2r(length - startPoint);

	if (ifPrintLog) {
		logger <<Priority::DEBUG <<"r=" <<r <<"\tr1=" <<r1 <<"\tr2=" <<r2 <<"\td1="
				<<startPoint <<"\td2=" <<(length - startPoint) <<"\td="
				<<length;
	}

	///计算基因型条件概率分布全数据，相当于生成书107页表
	vector<double> rMatrix(GENE_CP_SIZE);
	vector<string> geneMatrix(GENE_CP_SIZE);
	calcGeneCP(geneMatrix, rMatrix, f, m, q, r1, r2, r, 1, true, ifPrintLog);

	///统计各类型分布比率
	vector<vector<double> > gp = vector<vector<double> >(sampleSize, vector<double>(3));
	calculateGP(gp, mk, q, sampleSize, currentTrait, currentTrait + 1, geneMatrix, rMatrix, ifPrintLog);

	double u1 = 0, u2 = 0, u3 = 0, s1 = 0;
	/// EM算法
	EMCalculate(gp, sampleSize, expData, u0, s0, &u1, &u2, &u3, &s1, ifPrintLog);

	/// 计算LOD
	return calculateLOD(gp, sampleSize, expData, u0, s0, u1, u2, u3, s1, ifPrintLog);
}

/**
 * 打印帮助信息
 * @param fileName 文件名
 */
void printHelp(const string fileName) {
	cout << "Usage:\t" << fileName << " [ConfigFileName]" << endl;
}

/**
 * 运行一次QTL区间分析
 * @param args 参数个数
 * @param argv 参数内容
 * @return 最大LOD值
 */
double IMQTLRun(int args, char* argv[]) {
	int iSampleSize;
	int iTraitNumber;
	double step;

	string sExpressDataFile;
	string sChildrenGeneDataFile;
	string sFParentsGeneDataFile;
	string sMParentsGeneDataFile;
	string sTraitIntervalFile;

	bool ifPrintInitDataReport;
	bool ifPrintCalReport;
	bool ifPrintDetailedFinalReport;

	int LODThresholdTimer;

	if (args == 1) {
		/// 若参数个数为1，为手动模式运行，从操作台读取各个参数。
		inputValueFromKeyboard("the Number of Samples", &iSampleSize);
		inputValueFromKeyboard("the Number of Traits", &iTraitNumber);
		inputValueFromKeyboard("the length of step", &step);
		inputValueFromKeyboard("the Filename of Express Data", &sExpressDataFile);
		inputValueFromKeyboard("the Filename of Gene Data for Children", &sChildrenGeneDataFile);
		inputValueFromKeyboard("the Filename of Gene Data for Parent(F)", &sFParentsGeneDataFile);
		inputValueFromKeyboard("the Filename of Gene Data for Parent(M)", &sMParentsGeneDataFile);
		inputValueFromKeyboard("the Filename of Trait Interval Data", &sTraitIntervalFile);

		ifPrintInitDataReport = inputBooleanFromKeyboard("Do you want detailed Initialization Report?");
		ifPrintCalReport = inputBooleanFromKeyboard("Do you want detailed Calculation Report?");
		ifPrintDetailedFinalReport = inputBooleanFromKeyboard("Do you want detailed Final Report?");
	} else if (args >= 2) {
		/// 若参数个数大于2，为自动模式，从配置文件读取参数。
		logger <<Priority::INFO << "Input File:" << argv[1];
		fstream fin;
		fin.open(argv[1], ios::in);

		fin >> iSampleSize >> iTraitNumber >> step;
		fin >> sTraitIntervalFile;
		fin >> sChildrenGeneDataFile;
		fin >> sFParentsGeneDataFile >>sMParentsGeneDataFile;
		fin >> sExpressDataFile;
		fin >> ifPrintInitDataReport >> ifPrintCalReport >> ifPrintDetailedFinalReport;

		fin >> LODThresholdTimer;

		fin.close();
	} else {
		/// 参数个数错误，打印帮助。
		printHelp(argv[0]);
	}

	/// 汇总显示各输入参数
	logger <<Priority::INFO << "Input Summary:\t" << iSampleSize << " samples of " << iTraitNumber <<"\n"
			<< " traits in following files:" <<"\n"
			<< "\tTrait Interval Data in:     " << sTraitIntervalFile <<"\n"
			<< "\tGene Data of Children in:   " << sChildrenGeneDataFile <<"\n"
			<< "\tGene Data of Parent(F) in:  " << sFParentsGeneDataFile <<"\n"
			<< "\tGene Data of Parent(M) in:  " << sMParentsGeneDataFile <<"\n"
			<< "\tExpress Data in:            " << sExpressDataFile <<"\n"
			<< "For LOD Threshold, " <<LODThresholdTimer <<" random data will run." <<"\n"
			<< "Log Detail: " <<"\n"
			<< "\tInitialization Report - " << (ifPrintInitDataReport ? "On" : "Off") <<"\n"
			<< "\tCalculation Report    - " << (ifPrintCalReport ? "On" : "Off") <<"\n"
			<< "\tDetailed Final Report - " << (ifPrintDetailedFinalReport ? "On" : "Off") <<"\n"
			<< "Calculation step:\t" << step;

	logger <<Priority::INFO <<"========= INIT  =========";
	/// 初始化并读取位点距离数据
	vector<double> traitInterval(iTraitNumber);
	initIntervalData(sTraitIntervalFile, iTraitNumber - 1, traitInterval,
			ifPrintInitDataReport);

	/// 初始化基因型数据并读取基因型数据文件
	vector<vector<string> > mk;
	mk = vector<vector<string> >(iTraitNumber, vector<string>(iSampleSize, "  "));
	initChildrenGeneData(sChildrenGeneDataFile, mk, ifPrintInitDataReport);

	/// 初始化亲本的基因型
	vector<string> fGene(iTraitNumber, "abab");
	initParentGeneData(sFParentsGeneDataFile, fGene, iTraitNumber, "Parent(F)", ifPrintInitDataReport);
	vector<string> mGene(iTraitNumber, "abab");
	initParentGeneData(sMParentsGeneDataFile, mGene, iTraitNumber, "Parent(M)", ifPrintInitDataReport);
	vector<string> qGene(iTraitNumber, "Qq");

	/// 读取表型数据文件，并计算表型数组的均值、方差
	vector<EXPData> expData(iSampleSize);
	EXPData u0 = 0, s0 = 0;
	initExpressData(sExpressDataFile, iSampleSize, expData, &u0, &s0,
			ifPrintInitDataReport);

	logger <<Priority::DEBUG <<"init done";

	/// 逐个位点计算LOD值，并统计最大值
	double maxLOD = 0;
	int maxLODTrait;
	int maxLODPosition;
	logger <<Priority::INFO <<"========= START =========";
	for (int currentTrait = 0; currentTrait < iTraitNumber - 1;
			currentTrait++) {
		if (ifPrintDetailedFinalReport) {
			logger <<Priority::INFO << "---------" << (currentTrait + 1) << "---("
					<< traitInterval[currentTrait] << ")--------";
		}
		/// 按步长推进计算LOD值
		for (double startPoint = 0.0; startPoint < traitInterval[currentTrait];
				startPoint += step) {
			double LOD = intervalQTL(currentTrait, u0, s0,
					traitInterval[currentTrait], startPoint, mk, expData, iSampleSize, fGene[currentTrait], mGene[currentTrait], qGene[currentTrait], ifPrintCalReport);
			if (LOD > maxLOD) {
				maxLOD = LOD;
				maxLODTrait = currentTrait;
				maxLODPosition = startPoint;
			}
			if (ifPrintDetailedFinalReport) {
				logger <<Priority::INFO << "[" << startPoint << "~" << (startPoint + step) << "]:["
						<< traitInterval[currentTrait] << "] LOD=" << LOD;
			}
		}
		logger <<Priority::INFO <<"Trait " <<currentTrait+1 <<"~" <<currentTrait+2 <<" done\t("<<100*(currentTrait+1)/(iTraitNumber-1) <<"%)";
	}
	logger <<Priority::INFO <<"MAX LOD:\tTrait " <<maxLODTrait+1 << "[" << maxLODPosition << "~" << (maxLODPosition + step) << "] = " <<maxLOD;
	logger <<Priority::INFO <<"=========  END  =========";
	return maxLOD;
}

int main(int args, char* argv[]) {
	initLoggers();

	getTimeStamp();
	double re = IMQTLRun(args, argv);
	float t = getTimeStamp();

	logger <<Priority::INFO <<"Total time: " <<t <<" ms";
	logger <<Priority::INFO <<"Average time:" << t <<" ms";

	haltLoggers();
	return 0;
}
