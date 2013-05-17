#include "QTLUtils.h"

using namespace std;

/**
 * 获取系统时间的时间戳
 * @return 距离上一个时间戳的时间差，单位毫秒
 */
float getTimeStamp() {
	static clock_t t=clock();
	clock_t last = t;
	t = clock();
	return (float)(t-last)/ CLOCKS_PER_SEC * SECOND_TO_MS;
}

/**
 * 从键盘输入判断值，1为真，0为假
 * @param info 提示信息
 * @return 判断值
 */
int inputBooleanFromKeyboard(string info) {
	int re;
	cout <<info <<"(1.True, 0.False) ";
	cin >>re;
	cout <<(re ? "True" : "False") <<endl;
	return re;
}

/**
 * Haldane作图函数
 * @param d 遗传距离
 * @return 重组率
 */
double d2rHaldane(const double d) {
	return 0.5 * (1 - exp(-2 * d * 0.01));
}

/**
 * Kosambi作图函数
 * @param d 遗传距离
 * @return 重组率
 */
double d2rKosambi(const double d) {
	return 0.5 * tanh(2 * d * 0.01);
}

/**
 * 作图函数，默认使用Haldane作图函数
 * @param x 遗传距离
 * @return 重组率
 */
double d2r(const double x) {
	return d2rHaldane(x);
}

//vector<string> generateAllGeneType(const string inputGene) {
//	vector<string> re(3);
//	re[0] = inputGene[0]+inputGene[0];
//	re[1] = inputGene[0]+inputGene[1];
//	re[2] = inputGene[1]+inputGene[1];
//	return re;
//}

/**
 * 重排一个字符串内字符，如bdcae重排为abcde
 * @param input 需要重排的字符串
 * @return
 */
string reorderStr(const string input) {
	string re = input;
	sort(re.begin(), re.end());
	return re;
}
