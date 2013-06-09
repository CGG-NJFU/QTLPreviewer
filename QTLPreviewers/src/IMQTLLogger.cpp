/*
 * IMQTLLogger.cpp
 *
 *  Created on: 2013-5-16
 *      Author: yiqingxu
 */

#include "IMQTLLogger.h"

Category& logger = log4cpp::Category::getInstance("rootCategory");

/**
 * 初始化日志系统
 * @return 返回是否初始化成功
 */
bool initLoggers() {
	try {
		PropertyConfigurator::configure("config/log4cpp.properties");
	} catch (ConfigureFailure& f) {
		cerr << "log4cpp configuration failed:" <<endl
				<< f.what() <<endl;
		return false;
	}

	logger << Priority::DEBUG <<"log4cpp is ready";
	return true;
}

/**
 * 关闭日志系统
 */
void haltLoggers() {
	logger <<Priority::DEBUG <<"log4cpp is shutting down";
	Category::shutdown();
}

/**
 * 打印样本的表现型数据
 * @param data 样本表现型数组
 */
void printExpressData(const vector<EXPData> data) {
	printVectorData(data, "ExpressData", Priority::INFO);
}

/**
 * 打印遗传图谱的遗传区间距离数据
 * @param data 区间距离数据
 */
void printIntervalData(const vector<double> data) {
	printVectorData(data, "TraitIntervalData", Priority::INFO);
}

/**
 * 打印样本的基因型数据
 * @param data 样本基因型数据
 */
void printChildrenGeneData(const vector<vector<string> > data) {
	logger <<Priority::INFO <<"=======Children MarkData======";
	stringstream sstr;
	for (unsigned int i=0; i<data.size(); i++) {
		sstr <<endl;
		sstr <<"[" <<i <<"(" <<data[i].size() <<")" <<"]" <<"\t";
		for (unsigned int j=0; j<data[i].size(); j++ ) {
			sstr <<data[i][j] <<" ";
		}
	}
	logger <<Priority::INFO <<sstr.str();
}

/**
 * 打印亲本的基因型数据
 * @param data 亲本基因型数据
 * @param parentString 亲本名称
 */
void printParentGeneData(const vector<string> data, const string parentString) {
	printVectorData(data, parentString+" MarkData", Priority::INFO);
}

/**
 * 打印基因型条件概率分布全数据
 * @param geneMatrix 基因型
 * @param rMatrix 条件概率
 */
void printGeneCP(const vector<string>& geneMatrix, const vector<double> rMatrix) {
	logger <<Priority::DEBUG <<"ffmm:fm ->\tffmm:fm\tConditional Probability";
	for (unsigned int i=0; i<geneMatrix.size(); i++) {
		bitset<4> bs_up(i/4);  ///前四位为亲代基因型
		bitset<2> bs_down(i%4); ///后二位为子代基因型
		logger <<Priority::DEBUG <<bs_up <<":" <<bs_down <<" -> \t"
				<<geneMatrix[i].substr(0,4) <<":" << geneMatrix[i].substr(4,2) <<"\t"
				<<rMatrix[i];
	}
}

/**
 * 打印LOD判定结果
 * @param ifOK 判定条件
 * @param methodName 方法名称
 * @param timer 计数器
 * @param refLOD 参考LOD值
 * @param maxLOD LOD最大值
 */
void printLODJudgement(const bool ifOK, const string methodName, const unsigned int timer,
		const double refLOD, const double maxLOD) {
	if (ifOK) { ///若通过特殊性检查，则QTL识别成功
		logger << Priority::NOTICE << "Max of " << timer - 1
				<< " " <<methodName <<" re-mix data is " << refLOD
				<< " " <<(refLOD>maxLOD?">":"<") <<" "
				<< fabs(refLOD / maxLOD-1) * 100 << "%";
		logger << Priority::WARN << "QTL DETECTED: LOD:" << maxLOD;
		logger << Priority::NOTICE << "Please check full log for QTL detail "
				<< maxLOD;
	} else { ///否则，该位点QTL不存在
		logger << Priority::WARN
				<< "QTL detection need more verification, please check.";
	}
}
