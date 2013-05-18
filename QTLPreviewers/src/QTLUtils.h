/*
 * QTLUtils.h
 *
 *  Created on: 2013-5-15
 *      Author: yiqingxu
 */

#ifndef QTLUTILS_H_
#define QTLUTILS_H_

#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <Iterator>
#include <bitset>
#include <math.h>

// using log4cpp
#include <log4cpp/Category.hh>
#include <log4cpp/Appender.hh>
#include <log4cpp/FileAppender.hh>
#include <log4cpp/OstreamAppender.hh>
#include <log4cpp/Layout.hh>
#include <log4cpp/BasicLayout.hh>
#include <log4cpp/Priority.hh>
#include <log4cpp/PropertyConfigurator.hh>

#define SECOND_TO_MS 1000

typedef double EXPData; //表现型数据类型

using namespace std;
using namespace log4cpp;

extern Category& logger;

//---------------------------------
//----------- 字符串函数 -----------
//---------------------------------
string reorderStr(const string input);

/**
 * 将基本类型转换问字符串类型
 * @param value 基本类型值
 * @return 对应的字符串
 */
template <typename T>
string to_string(T const& value) {
    stringstream sstr;
    sstr << value;
    return sstr.str();
}

//---------------------------------
//----------- 遗传学函数 -----------
//---------------------------------
double d2r(const double x);
double d2rKosambi(const double d);
double d2rHaldane(const double d);
//vector<string> generateAllGeneType(const string inputGene);

//---------------------------------
//------------- 时间 --------------
//---------------------------------
float getTimeStamp();

//------------------------------------
//------------- 数学计算 --------------
//------------------------------------
/**
 * 获取向量中一组的算术平均值
 * @param begin 起始数据
 * @param end   结尾数据
 * @return 算数平均值
 */
template<class Iterator>
double getAverage(const Iterator begin, const Iterator end) {
	Iterator it;
	double sum = 0;
	for ( it = begin; it != end; it++) {
		sum += *it;
	}
	return sum/double(end-begin);
}

/**
 * 获取方差
 * @param begin 开始位置的游标
 * @param end 结束位置的游标
 * @param average 平均值，默认为0，视为自动计算
 * @return 数组数据的方差
 */
template <class Iterator>
double getVariance(const Iterator begin, const Iterator end, const double average) {
	double a = 0;
	if (0 == average) {
		a = getAverage(begin, end);
	}

	Iterator it;
	double sum = 0;
	int size = 0;
	for ( it = begin; it != end; it++, size++) {
		sum += (*it - average) * (*it - average);
	}

	return sum / (double)size;
}

//-----------------------------------
//------------- 键盘输入--------------
//-----------------------------------

/**
 * 从键盘输入变量的值
 * @param variableName 变量名
 * @param variableValue 变量值
 */
template <typename T>
void inputValueFromKeyboard(string variableName, T* variableValue) {
	cout <<"Please input " <<variableName <<":";
	cin >>*variableValue;
	cout <<variableName <<":" <<*variableValue;
	cout <<endl;
}

int inputBooleanFromKeyboard(string info=NULL);

//-----------------------------------
//------------- 屏幕输出--------------
//-----------------------------------
template <typename T>
void printVectorData(const vector<T> data, const string title = NULL) {
	if ( !title.empty() ) {
		logger <<Priority::INFO << "=======" <<title <<"=======";
	}
	for (unsigned int i = 0; i < data.size(); i++) {
		logger <<Priority::INFO << "[" << i << "]\t" << data[i];
	}
}

//-----------------------------------
//------------- 文件输入--------------
//-----------------------------------
/**
 * 从文件读入向量
 * @param filename 文件名
 * @param data 向量数据
 * @param size 向量大小，默认0为自动根据
 * @return 实际读入向量的大小
 */
template <typename T>
int readFile2Vector(string filename, vector<T>& data, const int size = 0) {
	fstream fin;
	fin.open(filename.data(), ios::in);

	if (size != 0) {
		data = vector<T>(size);
		for (int i=0; i<size; i++) {
			fin >> data[i];
		}
	} else {
		data = vector<T>();
		while (!fin.eof()) {
			T value;
			fin >> value;
			data.push_back(value);
		}
	}

	fin.close();
	return data.size();
}

/**
 * 从文件读入矩阵
 * @param filename 文件名
 * @param m 矩阵
 * @return
 */
template <typename T>
void readFile2Matrix(const string filename, vector< vector<T> >& m) {
	fstream fin;
	fin.open(filename.data(), ios::in);

	for (unsigned int i=0; i<m.size(); i++) {
		for (unsigned int j=0; j<m[i].size(); j++) {
			T in;
			fin >> in;
			m[i][j] = in;
		}
	}
	fin.close();
}

#endif /* QTLUTILS_H_ */
