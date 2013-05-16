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
#include <math.h>

#define SECOND_TO_MS 1000

typedef double EXPData; //表现型数据类型

using namespace std;

//---------------------------------
//----------- 遗传学函数 -----------
//---------------------------------
double d2r(double x);
double d2rKosambi(double d);
double d2rHaldane(double d);

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
 * @param data 浮点数组
 * @param size 数组尺寸
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
 * @param variableNamen 变量名
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
		cout << "=======" <<title <<"=======" <<endl;
	}
	for (unsigned int i = 0; i < data.size(); i++) {
		cout << "[" << i << "]\t" << data[i] <<"\n";
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

#endif /* QTLUTILS_H_ */
