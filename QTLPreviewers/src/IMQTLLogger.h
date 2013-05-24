/*
 * IMQTLLogger.h
 *
 *  Created on: 2013-5-16
 *      Author: yiqingxu
 */

#ifndef IMQTLLOGGER_H_
#define IMQTLLOGGER_H_

#include "QTLUtils.h"

using namespace std;
using namespace log4cpp;

void printExpressData(const vector<EXPData> data);
void printIntervalData(const vector<double> data);
void printChildrenGeneData(const vector<vector<string> > data);
void printParentGeneData(const vector<string> data, const string parentString);
void printGeneCP(const vector<string>& geneMatrix, const vector<double> rMatrix);
void printLODJudgement(const bool ifOK, const string methodName, const unsigned int timer,
		const double refLOD, const double maxLOD);

extern Category& logger;
bool initLoggers();
void haltLoggers();

#endif /* IMQTLLOGGER_H_ */
