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

void printExpressData(const vector<EXPData> data, const EXPData u0, const EXPData s0);
void printIntervalData(const vector<double> data);
void printGeneData(const vector<vector<string> > data);
void printGeneCP(const vector<string>& geneMatrix, const vector<double> rMatrix);

#endif /* IMQTLLOGGER_H_ */
