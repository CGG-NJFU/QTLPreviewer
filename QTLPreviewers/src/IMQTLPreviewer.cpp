//============================================================================
// Name        : IMQTLPreviewer.cpp
// Author      : Ethan
// Version     : Master
// Copyright   : Ethan @2012
// 运行配置文件见 dara/para.in
// 日志配置文件见 config/log4cpp.properties
// 主要编译选项见 src/IMQTLPreviewer.h
// Description : Previewer for Interval Mapping of QTL Analysis
//============================================================================

#include "IMQTLPreviewer.h"

extern Category& logger;

using namespace std;
using namespace log4cpp;

//缓存变量
vector<vector<string> > cachedGeneMatrix;
vector<vector<double> > cachedRMatrix;
vector<string> cachedFind;
vector<double> cachedResult;

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
 */
void initExpressData(const string fileName, const int sampleNumber, vector<EXPData>& expData) {
	logger <<Priority::DEBUG <<"Express data init start";
	int realsize = readFile2Vector(fileName, expData);
	if ( sampleNumber != realsize ) {
		logger <<Priority::WARN <<"Sample size (" <<realsize <<") might be wrong, please check.";
	}
	printExpressData(expData);
}

/**
 * 初始化样本的基因型数据
 * @param fileName 样本数据文件的文件名
 * @param data 数据数组
 */
void initChildrenGeneData(const string fileName, vector<vector<string> >& data) {
	logger <<Priority::DEBUG <<"Children Gene data init start";
	readFile2Matrix(fileName, data);
	printChildrenGeneData(data);
}

/**
 * 初始化亲本的基因型数据
 * @param fileName 亲本数据文件的文件名
 * @param geneData 数据数组
 * @param traitNumber 位点数量
 * @param parentName 亲本名称
 */
void initParentGeneData(const string fileName, vector<string>& geneData, const int traitNumber, const string parentName) {
	logger <<Priority::DEBUG <<"Parent Gene data init start";
	int realsize = readFile2Vector(fileName, geneData);
	if ( traitNumber != realsize ) {
		logger <<Priority::WARN <<"Trait number (" <<realsize <<") might be wrong, please check.";
	}

	printParentGeneData(geneData, parentName);
}

/**
 * 初始化遗传图谱的遗传区间距离数据
 * @param fileName 区间距离数据文件的文件名
 * @param intervalNumber 区间个数
 * @param data 区间距离数据
 */
void initIntervalData(const string fileName, const int intervalNumber, vector<double>& data) {
	logger <<Priority::DEBUG <<"Interval data init start";
	int realsize = readFile2Vector(fileName, data);
	if ( intervalNumber != realsize ) {
		logger <<Priority::WARN <<"The number of intervals (" <<realsize <<") might be wrong, please check.";
	}
	printIntervalData(data);
}

/**
 * 查找基因型的条件概率
 * @param geneMatrix 基因型
 * @param rMatrix 基因型概率
 * @param find 查找的匹配基因型
 * @return 基因型的条件概率
 */
double findGeneCP(const vector<string> geneMatrix, const vector<double> rMatrix, const string find, const bool ifUseCache=IF_USE_CACHE) {
	logger <<Priority::NOTSET <<"start to find Gene Conditional Possibility of " <<find
			<<", " <<(ifUseCache?"":"not ") << "using cache";

	bool ifInCacheFlag = true;

	///若启用缓存机制
	if (ifUseCache) {
		unsigned int i;
		///在缓存中查找
		for (i=0; i<cachedFind.size(); i++) {
			if ( cachedGeneMatrix[i]==geneMatrix &&
				cachedRMatrix[i]==rMatrix &&
				cachedFind[i]==find ){
				logger <<Priority::DEBUG <<"match cache";

				ifInCacheFlag = true;
				return cachedResult[i];
			}
		}

		///没有命中缓存，修改标记
		if (i==cachedFind.size()) {
			ifInCacheFlag = false;
		}
	}

	double same=0; //该基因型概率之和
	double all=0; //总概率之和，相当书上表格横行之和
	for (unsigned int i=0; i<geneMatrix.size(); i++) {
		if ( geneMatrix[i].substr(0,GENE_CP_CUT) == find.substr(0,GENE_CP_CUT) ) {
			logger <<Priority::NOTSET <<geneMatrix[i] <<" added:\t" <<rMatrix[i];
			all += rMatrix[i];

			if (geneMatrix[i]==find) {
				same += rMatrix[i];
			}
		}
	}

	double re;
	if ( all == 0 ) re = 0; else re = same/all;

	logger <<Priority::NOTSET <<"CP: " <<re <<" (" <<same <<"/" <<all <<")";

	if (ifUseCache) {
		if ( !ifInCacheFlag ) {
			logger <<Priority::DEBUG <<"not match cache, adding";
			///如果使用缓存，且结果不在缓存序列中，将结果加入缓存序列
			cachedGeneMatrix.push_back(geneMatrix);
			cachedRMatrix.push_back(rMatrix);
			cachedFind.push_back(find);
			cachedResult.push_back(re);
		}
	}

	return re;
}

/**
 * 根据基因型数据统计对应的实际分布概率
 * @param gp 实际分布概率
 * @param mk 基因型数据数组
 * @param qtlGene QTL基因的字符串
 * @param sampleIndex1 样本编号1
 * @param sampleIndex2 样本编号2
 * @param geneMatrix 基因型
 * @param rMatrix 条件型概率
 */
void calcGP(vector<vector<double> >& gp, const vector<vector<string> >& mk, const string qtlGene, const int sampleIndex1,
		const int sampleIndex2, const vector<string>& geneMatrix, const vector<double>& rMatrix) {
	///生成QTL基因型
	vector<string> qtl = generateAllGeneType(qtlGene);

	unsigned int sampleSize = mk[0].size();
	if (sampleSize != gp.size()) {
		logger <<Priority::WARN <<"DATASET size error might be in calculateGP function";
	}

	logger <<Priority::DEBUG <<"start to calculate Gene Possibility for all children.";

	for (unsigned int i = 0; i < sampleSize; i++) {
		///获取并拼接样本的基因型
		string gene = mk[sampleIndex1][i] + mk[sampleIndex2][i];

		///计算实际分布概率
		for (int j=0; j<3; j++) {
			gp[i][j] = findGeneCP(geneMatrix, rMatrix, gene+qtl[j], IF_USE_CACHE);
		}
	}

	///若开启缓存，清除缓存
	cachedGeneMatrix.clear();
	cachedRMatrix.clear();
	cachedFind.clear();
	cachedResult.clear();

	logger <<Priority::DEBUG << "=======GPData======";
	for (unsigned int c = 0; c < sampleSize; c++) {
		logger <<Priority::DEBUG << c << "\t" <<gp[c][0] <<"\t" <<gp[c][1] <<"\t"<<gp[c][2];
	}
}

/**
 * 利用LOD计分法计算LOD值
 * @param gp 实际分布概率
 * @param expData 表现型数据
 * @param s0 表现型的统计量sigima平方
 * @param s1
 * @param u0 表现型的统计量mu
 * @param u1
 * @param u2
 * @param u3
 * @return LOD计分
 */
double calcLOD(const vector<vector<double> > gp, const vector<EXPData> expData,
		const EXPData u0, const EXPData s0,
		const EXPData u1, const EXPData u2, const EXPData u3, const EXPData s1 ) {
	double sum1 = 0, sum2 = 0;

	unsigned int sampleSize = gp.size();
	if (sampleSize != expData.size()) {
		logger <<Priority::WARN <<"DATASET size error might be in calculateLOD function";
	}

	logger <<Priority::DEBUG <<"start to calculate LOD";
	logger <<Priority::DEBUG <<"u0=" <<u0 <<"\ts0=" <<s0;
	for (unsigned int i = 0; i < sampleSize; i++) {
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

	logger <<Priority::INFO <<"sum1=" <<sum1 <<"\tsum2=" <<sum2 <<"\tLOD=" <<lod;

	return lod;
}

/**
 * 利用EM算法迭代计算最大似然估计
 * @param expData 表现型数据
 * @param gp 实际分布概率
 * @param u0 样本的统计量mu
 * @param s0 样本的统计量sigma平方
 * @param u1
 * @param u2
 * @param u3
 * @param s1
 */
void EMCalc(const vector<vector<double> >& gp, const vector<EXPData> expData,
		const EXPData u0, const EXPData s0, EXPData* u1, EXPData* u2, EXPData* u3, EXPData* s1) {
	unsigned int sampleSize = expData.size();
	if (sampleSize != gp.size()) {
		logger <<Priority::WARN <<"DATASET size error might be in EMCalculate function";
	}

	/// 计算u1, u2, u3, sigma2
	EXPData u10 = u0;
	EXPData u20 = u0;
	EXPData u30 = u0;
	EXPData s10 = s0;	//*s1+1;

	/// 提取收敛精度要求
	int step = 1; //日志输出分步初试值为1

	logger <<Priority::DEBUG <<"=======EMSteps======";
	logger <<Priority::DEBUG <<"[step]\tu1\tu2\tu3\ts1";

	while (fabs(u10 - *u1) > ACCURACY_REQ || fabs(u20 - *u2) > ACCURACY_REQ
			|| fabs(u30 - *u3) > ACCURACY_REQ || fabs(s10 - *s1) > ACCURACY_REQ) {
		double py1 = 0, py10 = 0;
		double py2 = 0, py20 = 0;
		double py3 = 0, py30 = 0;
		double py4 = 0, py5 = 0, py6 = 0;

		unsigned int i;
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

		logger <<Priority::DEBUG << "[" << step++ << "]\t" << *u1 << "\t" << *u2 << "\t" << *u3
					<< "\t" << *s1;
	}	//end-while

	//输出终值
	*u1 = u10;
	*u2 = u20;
	*u3 = u30;
	*s1 = s10;
	logger <<Priority::INFO << "Result:\tu1=" <<*u1 <<"\tu2=" <<*u2 <<"\tu3=" <<*u3 <<"\ts1=" <<*s1;
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
void sortGene(vector<string>& geneMatrix) {
	logger <<Priority::DEBUG <<"Reordering Gene...";
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
 * @return 成功运行则返回1
 */
int calcGeneCP(vector<string>& geneMatrix, vector<double>& rMatrix,
		const string f, const string m, const string q,
		const double r1, const double r2, const double r=0, const double zeta=1,
		const bool ifReorder=true) {

	unsigned int mSize = geneMatrix.size();
	if (mSize != rMatrix.size()) {
		logger <<Priority::WARN <<"DATASET size error might be in calcGeneCP function";
	}

	logger <<Priority::DEBUG <<"start to calculate full Conditional Possibility for all Children";
	for (unsigned int i=0; i<mSize; i++) {
		bitset<4> bs_fm(i/GENE_CP_CUT); //前4位为亲本基因型
		bitset<2> bs_qtl(i%GENE_CP_CUT); //后2位为QTL基因型

		int fGene = bs_fm[3]*4 + bs_qtl[1]*2 + bs_fm[1]; //父本基因型
		int mGene = bs_fm[2]*4 + bs_qtl[0]*2 + bs_fm[0]; //母本基因型

		bitset<3> bs_up(fGene);
		bitset<3> bs_down(mGene);

		/// 计算并填充条件概率
		rMatrix[i] = calcCross(bs_up, r1, r2, r, zeta) * calcCross(bs_down, r1, r2, r, zeta);

		/// 计算并填充基因型数据
		geneMatrix[i] =
				(i&0x20?f.substr(1,1):f.substr(0,1)) +
				(i&0x10?f.substr(3,1):f.substr(2,1)) +
				(i&0x08?m.substr(1,1):m.substr(0,1)) +
				(i&0x04?m.substr(3,1):m.substr(2,1)) +
				(i&0x02?q.substr(1,1):q.substr(0,1)) +
				(i&0x01?q.substr(1,1):q.substr(0,1));
	}

	if (ifReorder) {
		sortGene(geneMatrix);
	}

	printGeneCP(geneMatrix, rMatrix);

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
 * @return 返回LOD值
 */
double intervalQTL(const int currentTrait, const EXPData u0, const EXPData s0, const double length,
		const double startPoint, const vector<vector<string> >& mk, const vector<EXPData> expData, const int sampleSize,
		const string f, const string m, const string q) {
	double r = d2r(length);
	double r1 = d2r(startPoint);
	double r2 = d2r(length - startPoint);
	logger <<Priority::DEBUG <<"start to interval QTL:";
	logger <<Priority::DEBUG <<"r=" <<r <<"\tr1=" <<r1 <<"\tr2=" <<r2 <<"\td1="
			<<startPoint <<"\td2=" <<(length - startPoint) <<"\td="
			<<length;

	///计算基因型条件概率分布全数据，相当于生成书107页表
	vector<double> rMatrix(GENE_CP_SIZE);
	vector<string> geneMatrix(GENE_CP_SIZE);
	calcGeneCP(geneMatrix, rMatrix, f, m, q, r1, r2, r, 1, true);

	///统计各类型分布比率
	vector<vector<double> > gp = vector<vector<double> >(sampleSize, vector<double>(3));
	calcGP(gp, mk, q, currentTrait, currentTrait + 1, geneMatrix, rMatrix);

	double u1 = 0, u2 = 0, u3 = 0, s1 = 0;
	/// EM算法
	EMCalc(gp, expData, u0, s0, &u1, &u2, &u3, &s1);

	/// 计算LOD
	return calcLOD(gp, expData, u0, s0, u1, u2, u3, s1);
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
 * @param mk 子代基因型数据
 * @param expData 子代表现型数据
 * @param traitInterval 位点间隔距离
 * @param fGene 父本基因型数据
 * @param mGene 母本基因型数据
 * @param qGene QTL基因型数据
 * @param step 步长
 * @param iRemixExpDataTimer 随机检测的计数器，0表示使用原始数据。
 * @return 最大LOD值
 */
double QTLRun(const vector<vector<string> >& mk, const vector<EXPData>& expData, const vector<double>& traitInterval,
		const vector<string>& fGene, const vector<string>& mGene, const vector<string>& qGene,
		const double step, const int iRemixExpDataTimer = 0) {
	int iTraitNumber = traitInterval.size();
	unsigned int iSampleSize = mk[0].size();
	if (iSampleSize != expData.size()) {
		logger <<Priority::WARN <<"DATASET size error might be in QTLRun function";
	}

	/// 初始化表型数据，计算算术平均值和方差
	EXPData u0 = getAverage( expData.begin(), expData.end() );
	EXPData s0 = getVariance( expData.begin(), expData.end(), u0);

	vector<EXPData> useExpData;

	if (iRemixExpDataTimer==0) { //重排模式下不输出详细日志
		logger <<Priority::INFO <<"u0=" <<u0 << "\ts0=" <<s0;
	}
	useExpData = expData;

	/// 逐个位点计算LOD值，并统计最大值
	double maxLOD = 0;
	int maxLODTrait;
	int maxLODPosition;

	for (int currentTrait = 0; currentTrait <= iTraitNumber - 1;
			currentTrait++) {
		if (iRemixExpDataTimer==0) { //重排模式下不输出详细日志
			logger <<Priority::INFO << "---------" << (currentTrait + 1) << "---("
					<< traitInterval[currentTrait] << ")--------";
		}

		/// 按步长推进计算LOD值
		for (double startPoint = 0.0; startPoint < traitInterval[currentTrait];
				startPoint += step) {
			double LOD = intervalQTL(currentTrait, u0, s0,
					traitInterval[currentTrait], startPoint, mk, useExpData, iSampleSize, fGene[currentTrait]+fGene[currentTrait+1], mGene[currentTrait]+mGene[currentTrait+1], qGene[currentTrait]);
			if (LOD > maxLOD) {
				maxLOD = LOD;
				maxLODTrait = currentTrait;
				maxLODPosition = startPoint;
			}

			if (iRemixExpDataTimer==0) { //重排模式下不输出详细日志
				logger <<Priority::NOTICE <<currentTrait+1 <<"\t[ " << startPoint << " ~ "
						<< ((startPoint+step)>traitInterval[currentTrait]?traitInterval[currentTrait]:(startPoint+step))
						<< " ]:[ " << traitInterval[currentTrait] << " ]\t" << LOD;
			}
		}
		if (iRemixExpDataTimer==0) { //重排模式下不输出详细日志
			logger <<Priority::INFO <<"Trait " <<currentTrait+1 <<"~" <<currentTrait+2 <<" done\t("<<100*(currentTrait+1)/(iTraitNumber-1) <<"%)";
		}
	}
	if (iRemixExpDataTimer==0) { //重排模式下不输出详细日志
		logger <<Priority::NOTICE <<"MAX LOD:\tTrait " <<maxLODTrait+1 << "[" << maxLODPosition << "~"
				<< (maxLODPosition+step>traitInterval[maxLODTrait]?traitInterval[maxLODTrait]:maxLODPosition+step)
				<< "] = " <<maxLOD;
	} else {
		logger <<Priority::NOTICE <<"[" <<iRemixExpDataTimer <<"] MAX LOD in Random Re-mix: = " <<maxLOD;
	}
	return maxLOD;
}

/**
 * 判定LOD参考值与假定LOD最大值的关系
 * @param refLOD LOD参考值
 * @param maxLOD 假定LOD最大值
 * @return 若参考值超过阈值，则返回True，否则返回false
 */
bool judgeLOD(const double refLOD, const double maxLOD) {
	/// 如果LOD参考值大于假定LOD最大值，打印警报日志
	if (refLOD >= maxLOD) {
		logger << Priority::WARN << "maxLOD (" << maxLOD
				<< ") is challenged by a re-mix data[" << refLOD << ">(" <<(refLOD/maxLOD-1)*100 <<"%)]";
		/// 如果LOD参考值大于假定LOD最大值的阈值，退出，判定LOD无效
		if (refLOD >= maxLOD * LOD_THRESHOLD)
			return true;
	}
	return false;
}

/**
 * 利用排列数向量混淆表型数据，计算LOD参考值
 * @param permutationOrder 排列数首序列标号
 * @param max_permutation 最大排列数
 * @param maxLOD 假定LOD最大值
 * @param mk 子代基因型数据
 * @param expData 子代表现型数据
 * @param traitInterval 位点间隔距离
 * @param fGene 父本基因型数据
 * @param mGene 母本基因型数据
 * @param qGene QTL基因型数据
 * @param dStep 步长
 * @return 总计算次数
 */
unsigned int LODPermutationRun(const int permutationOrder,
		const double maxLOD, const vector<vector<string> >& mk,
		const vector<EXPData>& expData, const vector<double>& traitInterval,
		const vector<string>& fGene, const vector<string>& mGene,
		const vector<string>& qGene, const double dStep, const unsigned int max_permutation = MAX_PERMUTATION) {
	unsigned int timer;

	///启动LOD特殊性检查
	logger <<Priority::INFO <<"========= REMIX =========";

	///生成重排使用的排列数向量
	int sampleSize = expData.size();
	vector<int> order = first_permutation_order(sampleSize, permutationOrder);
	vector<int> end = last_permutation_order(sampleSize, permutationOrder);

	double secLOD=0;
	for (timer=1; order!=end && timer<=max_permutation; timer++) {
		///按排列数向量重排表型数据
		vector<EXPData> useExpData = vector<EXPData>(expData);
		reorderTargetBy(useExpData, order);

		///计算重排的LOD参考值
		double permutationLOD = QTLRun(mk, useExpData, traitInterval, fGene, mGene, qGene, dStep, timer);
		if (permutationLOD > secLOD) {
			secLOD = permutationLOD;
		}
		///判定LOD参考值
		if ( judgeLOD(permutationLOD, maxLOD) ) break;
		next_permutation(order.begin(), order.end());

	}
	///打印判定结果
	printLODJudgement( order==end||timer>max_permutation, "permutation", timer, secLOD, maxLOD);
	return timer-1;
}

/**
 * 利用随机排列混淆表型数据，计算LOD参考值
 * @param iRandomLODTimer 随机次数
 * @param maxLOD 假定LOD最大值
 * @param mk 子代基因型数据
 * @param expData 子代表现型数据
 * @param traitInterval 位点间隔距离
 * @param fGene 父本基因型数据
 * @param mGene 母本基因型数据
 * @param qGene QTL基因型数据
 * @param dStep 步长
 * @return 总计算次数
 */
unsigned int LODRandomRun(const unsigned int iRandomLODTimer, const double maxLOD,
		const vector<vector<string> >& mk, const vector<EXPData>& expData,
		const vector<double>& traitInterval, const vector<string>& fGene,
		const vector<string>& mGene, const vector<string>& qGene,
		const double dStep) {
	double secLOD = 0;
	unsigned int timer=1;
	///启动LOD特殊性检查
	if ( iRandomLODTimer > 0 ) {
		logger <<Priority::INFO <<"========= REMIX =========";

		for (timer = 1; timer<=iRandomLODTimer; timer++) {
			///若启用随机重排模式，随机重排所有表型数据
			vector<EXPData> useExpData;
			useExpData = vector<EXPData>(expData);
			random_shuffle(useExpData.begin(), useExpData.end());

			///计算LOD参考值
			double randomLOD = QTLRun(mk, useExpData, traitInterval, fGene, mGene, qGene, dStep, timer);
			if ( randomLOD > secLOD ) {
				secLOD = randomLOD;
			}
			/// 判定LOD参考值
			if ( judgeLOD(randomLOD, maxLOD) ) break;
		}
		///打印判定结果
		printLODJudgement( timer>iRandomLODTimer, "random", timer, secLOD, maxLOD);
	} else {
		///若配置中检查次数为0，则视为跳过检查阶段
		logger <<Priority::WARN <<"Re-mix verification passed, QTL detection is in-completed.";
	}
	return timer-1;
}

int mainRun(int args, char* argv[]) {
	int iSampleSize;
	int iTraitNumber;
	double dStep;

	string sExpressDataFile;
	string sChildrenGeneDataFile;
	string sFParentsGeneDataFile;
	string sMParentsGeneDataFile;
	string sTraitIntervalFile;

	bool bIfUseRandomLOD;
	unsigned int iRandomLODTimer = 0;
	int iPermutationOrder;

	/// 初始化日志系统、计时器开始计时
	bool ifLogOK = initLoggers();
	getTimeStamp();

	if (args == 1) {
		/// 若参数个数为1，为手动模式运行，从操作台读取各个参数。
		inputValueFromKeyboard("the Number of Samples", &iSampleSize);
		inputValueFromKeyboard("the Number of Traits", &iTraitNumber);
		inputValueFromKeyboard("the length of step", &dStep);
		inputValueFromKeyboard("the Filename of Express Data", &sExpressDataFile);
		inputValueFromKeyboard("the Filename of Gene Data for Children", &sChildrenGeneDataFile);
		inputValueFromKeyboard("the Filename of Gene Data for Parent(F)", &sFParentsGeneDataFile);
		inputValueFromKeyboard("the Filename of Gene Data for Parent(M)", &sMParentsGeneDataFile);
		inputValueFromKeyboard("the Filename of Trait Interval Data", &sTraitIntervalFile);

		bIfUseRandomLOD = inputBooleanFromKeyboard("Use Random LOD detection?");
		if (bIfUseRandomLOD) {
			inputValueFromKeyboard("the times of Random LOD Calculation", &iRandomLODTimer);
		} else {
			inputValueFromKeyboard("the Filename of Permutation Order", &iPermutationOrder);
		}
	} else if (args == 2) {
		/// 若参数个数大于2，为自动模式，从配置文件读取参数。
		logger <<Priority::INFO << "Input File:" << argv[1];
		fstream fin;
		fin.open(argv[1], ios::in);

		fin >> iSampleSize >> iTraitNumber >> dStep;
		fin >> sTraitIntervalFile;
		fin >> sChildrenGeneDataFile;
		fin >> sFParentsGeneDataFile >>sMParentsGeneDataFile;
		fin >> sExpressDataFile;

		fin >> bIfUseRandomLOD;
		if (bIfUseRandomLOD) {
			fin >> iRandomLODTimer;
		} else {
			fin >> iPermutationOrder;
		}

		fin.close();
	} else {
		/// 参数个数错误，打印帮助。
		printHelp(argv[0]);
		return 0;
	}

	/// 汇总显示各输入参数
	logger <<Priority::NOTICE << "\nInput Summary:\t" << iSampleSize << " samples of " << iTraitNumber
			<< " traits in following files:" <<"\n"
			<< "\tTrait Interval Data in:     " << sTraitIntervalFile <<"\n"
			<< "\tGene Data of Children in:   " << sChildrenGeneDataFile <<"\n"
			<< "\tGene Data of Parent(F) in:  " << sFParentsGeneDataFile <<"\n"
			<< "\tGene Data of Parent(M) in:  " << sMParentsGeneDataFile <<"\n"
			<< "\tExpress Data in:            " << sExpressDataFile <<"\n"
			<< "Calculation dStep:\t" << dStep;

	logger <<Priority::INFO <<"========= INIT  =========";
	/// 初始化并读取位点距离数据
	vector<double> traitInterval(iTraitNumber);
	initIntervalData(sTraitIntervalFile, iTraitNumber - 1, traitInterval);

	/// 初始化基因型数据并读取基因型数据文件
	vector<vector<string> > mk;
	mk = vector<vector<string> >(iTraitNumber, vector<string>(iSampleSize, "  "));
	initChildrenGeneData(sChildrenGeneDataFile, mk);

	/// 初始化亲本的基因型
	vector<string> fGene(iTraitNumber, "ab");
	initParentGeneData(sFParentsGeneDataFile, fGene, iTraitNumber, "Parent(F)");
	vector<string> mGene(iTraitNumber, "ab");
	initParentGeneData(sMParentsGeneDataFile, mGene, iTraitNumber, "Parent(M)");
	vector<string> qGene(iTraitNumber, "Qq");

	/// 读取表型数据文件，并计算表型数组的均值、方差
	vector<EXPData> expData(iSampleSize);
	initExpressData(sExpressDataFile, iSampleSize, expData);

	logger <<Priority::DEBUG <<"init done";
	logger <<Priority::INFO <<"========= START =========";
	/// 计算LOD最大值
	double maxLOD = QTLRun(mk, expData, traitInterval, fGene, mGene, qGene, dStep);

	/// 验证LOD最大值的有效性
	logger <<Priority::INFO <<"==== LOD Calculation ====";
	unsigned int timer=1;
	if (bIfUseRandomLOD) {
		/// 方法一：利用随机排列混淆表型数据，计算LOD参考值
		logger <<Priority::NOTICE <<"\nFor LOD Detection: " <<iRandomLODTimer <<" random data will run.";
		timer += LODRandomRun(iRandomLODTimer, maxLOD, mk, expData, traitInterval, fGene, mGene, qGene, dStep);
	} else {
		/// 方法二：利用排列数向量混淆表型数据，计算LOD参考值
		logger <<Priority::NOTICE <<"\nFor LOD Detection: OrderGroup \"" <<iPermutationOrder <<"/"<<iSampleSize <<"\" will run.";
		timer += LODPermutationRun(iPermutationOrder, maxLOD, mk, expData, traitInterval, fGene, mGene, qGene, dStep);
	}

	/// 若日志系统故障，只通过标准流输出最终结果，并提示检查日志系统
	if ( !ifLogOK ) {
		cout <<"maxLOD=" <<maxLOD <<endl
				<<"For more detail, please install and configure log4cpp in your system.";
	}
	logger <<Priority::INFO <<"=========  END  =========";

	/// 完成计时
	float t = getTimeStamp();
	logger <<Priority::NOTICE <<"Total time: " <<t <<" ms";
	logger <<Priority::NOTICE <<"Average time: " <<t/timer <<" ms";

	/// 退出日志系统
	haltLoggers();

	return 0;
}

int main(int args, char* argv[]) {
	return mainRun(args, argv);
}
