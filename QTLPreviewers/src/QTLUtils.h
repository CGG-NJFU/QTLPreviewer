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
 * 获取平均值
 * @param data 浮点数组
 * @param size 数组尺寸
 * @return 数组数据的平均值
 */
template <typename T>
T getAverage(T* data, int size) {
	int i;
	T sum = 0;
	for (i = 0; i < size; i++)
		sum += data[i];
	return sum / (double)size;
}

/**
 * 获取方差
 * @param data 浮点数组
 * @param size 数组尺寸
 * @param average 平均值，默认为0，视为自动计算
 * @return 数组数据的方差
 */
template <typename T>
T getVariance(T* data, int size, T average) {
	int i;
	T sum = 0;
	if (0 == average) {
		average = getAverage(data, size);
	}
	for (i = 0; i < size; i++) {
		sum += (data[i] - average) * (data[i] - average);
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

#endif /* QTLUTILS_H_ */
