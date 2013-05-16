/*
 * IMQTLLogger.cpp
 *
 *  Created on: 2013-5-16
 *      Author: yiqingxu
 */

#include "IMQTLLogger.h"

/**
 * 打印样本的表现型数据
 * @param data 样本表现型数组
 * @param u0 统计量mu
 * @param s0 统计量sigma平方
 */
void printExpressData(const vector<EXPData> data, const EXPData u0, const EXPData s0) {
	printVectorData(data, "ExpressData");
	cout << "u0=" << u0 << endl << "s0=" << s0 << endl;
}

/**
 * 打印遗传图谱的遗传区间距离数据
 * @param data 区间距离数据
 */
void printIntervalData(const vector<double> data) {
	printVectorData(data, "TraitIntervalData");
}
