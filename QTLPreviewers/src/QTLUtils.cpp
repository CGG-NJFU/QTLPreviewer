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
bool inputBooleanFromKeyboard(string info) {
	bool re;
	cout <<info <<" (1.True, Any Other.False) ";
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

/**
 * 生成全基因型排列，如Qq=》{QQ, Qq, qq}
 * @param inputGene 基因型字母
 * @return 基因型字母全排列
 */
vector<string> generateAllGeneType(const string inputGene) {
	vector<string> re(3);
	re[0] = inputGene.substr(0,1)+inputGene.substr(0,1); //"QQ"
	re[1] = inputGene.substr(0,1)+inputGene.substr(1,1); //"Qq"
	re[2] = inputGene.substr(1,1)+inputGene.substr(1,1); //"qq"
	return re;
}

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

/**
 * 生成序列顺序向量
 * @param length 序列长度
 * @param first 首序列标号，下标范围[0~length-1]
 * @return 生成的序列顺序向量
 */
vector<int> first_permutation_order(const unsigned int length, const unsigned int first) {
	vector<int> order(length);
	order[0] = (first-1)%length;
	for(unsigned int i=1;i<length;i++) {
		if (i<=(first-1)%length) order[i+1]=i;
		else order[i]=i;
	}
	return order;
}

vector<int> last_permutation_order(const unsigned int length, const unsigned int first) {
	vector<int> order;
	order = first_permutation_order(length, first%length+1);
	prev_permutation(order.begin(), order.end());
	return order;
}
