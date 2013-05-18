/*
 * IMQTLLogger.cpp
 *
 *  Created on: 2013-5-16
 *      Author: yiqingxu
 */

#include "IMQTLLogger.h"

Category& logger = log4cpp::Category::getInstance("rootAppender");

void initLoggers() {
	try {
		PropertyConfigurator::configure("config/log4cpp.properties");
	} catch (ConfigureFailure& f) {
		cerr << "log4cpp configuration failed:" <<endl
				<< f.what() << std::endl;
	}

	logger << Priority::DEBUG <<"LOG4CPP is ready";
}

void haltLoggers() {
	Category::shutdown();
}

/**
 * 打印样本的表现型数据
 * @param data 样本表现型数组
 * @param u0 统计量mu
 * @param s0 统计量sigma平方
 */
void printExpressData(const vector<EXPData> data, const EXPData u0, const EXPData s0) {
	printVectorData(data, "ExpressData");
	logger <<Priority::DEBUG <<"u0=" <<u0 << "\ts0=" <<s0;
}

/**
 * 打印遗传图谱的遗传区间距离数据
 * @param data 区间距离数据
 */
void printIntervalData(const vector<double> data) {
	printVectorData(data, "TraitIntervalData");
}

/**
 * 打印样本的基因型数据
 * @param data 样本基因型数据
 */
void printChildrenGeneData(const vector<vector<string> > data) {
	logger <<Priority::DEBUG <<"=======Children MarkData======";
	stringstream sstr;
	for (unsigned int i=0; i<data.size(); i++) {
		sstr <<endl;
		sstr <<"[" <<i <<"(" <<data[i].size() <<")" <<"]" <<"\t";
		for (unsigned int j=0; j<data[i].size(); j++ ) {
			sstr <<data[i][j] <<" ";
		}
	}
	logger <<Priority::DEBUG <<sstr.str();
}

/**
 * 打印亲本的基因型数据
 * @param data 亲本基因型数据
 * @param parentString 亲本名称
 */
void printParentGeneData(const vector<string> data, const string parentString) {
	printVectorData(data, parentString+" MarkData");
}

/**
 * 打印基因型条件概率分布全数据
 * @param geneMatrix 基因型
 * @param rMatrix 条件概率
 */
void printGeneCP(const vector<string>& geneMatrix, const vector<double> rMatrix) {
	logger <<Priority::DEBUG <<"ffmm:fm ->\tffmm:fm\tConditional Probability\n";
	for (unsigned int i=0; i<geneMatrix.size(); i++) {
		bitset<4> bs_up(i/4);  ///前四位为亲代基因型
		bitset<2> bs_down(i%4); ///后二位为子代基因型
		logger <<Priority::DEBUG <<bs_up <<":" <<bs_down <<" -> \t"
				<<geneMatrix[i].substr(0,4) <<":" << geneMatrix[i].substr(4,2) <<"\t"
				<<rMatrix[i];
	}
}
