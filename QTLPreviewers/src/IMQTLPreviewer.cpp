//============================================================================
// Name        : IMQTLPreviewer.cpp
// Author      : Ethan
// Version     : Master
// Copyright   : Ethan @2012
//		iSampleSize = 189;
//		iTraitNumber = 9;
//		step=1.0;
//
//		sExpressDataFile = "data/express-66.txt";
//		sGeneDataFile = "data/gene-66.txt";
//		sTraitIntervalFile = "data/traitInterval-66.txt";
//
//		ifPrintInitDataReport = false;
//		ifPrintCalReport = false;
//		ifPrintFinalReport = true;
// Description : Previewer for Interval Mapping of QTL Analysis
//============================================================================

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define ACCURACY_REQ 0.000001
#define VERBOSE_MODE 1

//概率密度函数
double f(double x, double u, double s2) {
	return ((1 / sqrt(2 * M_PI * s2)) * exp(-0.5 * (x - u) * (x - u) / s2));
}

//获取平均值
double getAverage(double* data, int size) {
	int i;
	double sum = 0;
	for (i = 0; i < size; i++)
		sum += data[i];
	return sum / size;
}

//获取方差
double getVariance(double* data, int size, double average = 0) {
	int i;
	double sum = 0;
	if (0 == average) {
		average = getAverage(data, size);
	}
	for (i = 0; i < size; i++)
		sum = sum + (data[i] - average) * (data[i] - average);
	return sum / size;
}

//打印表现型数据
void printExpressData(double *data, int size, double u0, double s0) {
	cout<<"=======ExpressData======="<<endl;
	for (int i = 0; i < size; i++) {
		cout<<"["<<i<<"]\t"<<data[i]<<endl;
	}
	cout <<"u0=" <<u0 <<endl <<"s0=" <<s0 <<endl;
}

//初始化表现型数据
void initExpressData(const string fileName, const int sampleNumber,
		double* data, double* u0, double* s0, int ifPrint=VERBOSE_MODE) {
	fstream fin;
	fin.open(fileName.data(),ios::in);
	for (int i=0; i<sampleNumber; i++) {
		fin >> data[i];
	}
	fin.close();

	*u0 = getAverage(data, sampleNumber);
	*s0 = getVariance(data, sampleNumber, *u0);
	if (ifPrint) {
		printExpressData(data, sampleNumber, *u0, *s0);
	}
}

//打印基因数据
void printGeneData(char** _data, int sampleNumber, int traitNumber) {
	char* data = (char*)_data;
	cout<<"=======MarkData======"<<endl;
	int i, j;
	for (i = 0; i < traitNumber; i++) {
		cout<<"["<<i<<"]\t";
		for (j = 0; j < sampleNumber; j++) {
			if (j%20==0) cout<<endl <<"[" <<j <<"]";
			cout<<*(data+i*sampleNumber+j);
		}
		cout<<endl;
	}
}

//初始化基因型数据
void initGeneData(string fileName, const int sampleNumber, const int traitNumber,
		char** _data, int ifPrint=VERBOSE_MODE) {
	char* data = (char*)_data;
	fstream fin;
	fin.open(fileName.data(),ios::in);

	for (int i=0; i<traitNumber; i++) {
		string line;
		getline(fin, line);

		for (int j=0; j<sampleNumber; j++) {
			*(data+i*sampleNumber+j) = line[j];
		}
	}
	fin.close();

	if(ifPrint) {
		printGeneData(_data, sampleNumber, traitNumber);
	}
}

//打印位点间隔数据
void printIntervalData(double* data, int intervalNumber) {
	cout<<"=======TraitIntervalData======="<<endl;
	for (int i = 0; i < intervalNumber; i++) {
		cout<<"["<<i<<"]\t"<<data[i]<<endl;
	}
}

//初始化位点间隔数据
void initIntervalData(string fileName, const int intervalNumber,
		double* data, int ifPrint=VERBOSE_MODE) {
	fstream fin;
	fin.open(fileName.data(),ios::in);
	for (int i=0; i<intervalNumber; i++) {
		fin >> data[i];
	}
	fin.close();

	if (ifPrint) {
		printIntervalData(data, intervalNumber);
	}
}

//计算Pij
void calculatePij(double p[9][3], double r, double r1, double r2, int ifPrint=VERBOSE_MODE) {
	//无干扰时F2群体 QTL 基因频率分布数组函数p 书107页//
	p[0][0] = (1 - r1) * (1 - r1) * (1 - r2) * (1 - r2) / ((1 - r) * (1 - r));
	p[0][1] = 2 * r1 * r2 * (1 - r1) * (1 - r2) / ((1 - r) * (1 - r));
	p[0][2] = r1 * r1 * r2 * r2 / ((1 - r) * (1 - r));
	p[1][0] = (1 - r1) * (1 - r1) * r2 * (1 - r2) / (r * (1 - r));
	p[1][1] = r1 * (1 - r1) * (1 - 2 * r2 + 2 * r2 * r2) / (r * (1 - r));
	p[1][2] = r1 * r1 * r2 * (1 - r2) / (r * (1 - r));
	p[2][0] = (1 - r1) * (1 - r1) * r2 * r2 / (r * r);
	p[2][1] = 2 * r1 * r2 * (1 - r1) * (1 - r2) / (r * r);
	p[2][2] = r1 * r1 * (1 - r2) * (1 - r2) / (r * r);
	p[3][0] = r1 * (1 - r1) * (1 - r2) * (1 - r2) / (r * (1 - r));
	p[3][1] = (1 - 2 * r1 + 2 * r1 * r1) * r2 * (1 - r2) / (r * (1 - r));
	p[3][2] = r1 * (1 - r1) * r2 * r2 / (r * (1 - r));
	p[4][0] = 2 * r1 * r2 * (1 - r1) * (1 - r2) / (1 - 2 * r + 2 * r * r);
	p[4][1] = (1 - 2 * r1 + 2 * r1 * r1) * (1 - 2 * r2 + 2 * r2 * r2)
			/ (1 - 2 * r + 2 * r * r);
	p[4][2] = p[4][0];
	p[5][0] = p[3][2];
	p[5][1] = p[3][1];
	p[5][2] = p[3][0];
	p[6][0] = p[2][2];
	p[6][1] = p[2][1];
	p[6][2] = p[2][0];
	p[7][0] = p[1][2];
	p[7][1] = p[1][1];
	p[7][2] = p[1][0];
	p[8][0] = p[0][2];
	p[8][1] = p[0][1];
	p[8][2] = p[0][0];

	if (ifPrint) {
		for (int i = 0; i < 9; i++) {
			double sum = 0;
			for (int j = 0; j < 3; j++) {
				sum += p[i][j];
				cout<<p[i][j]<<"\t";
			}
			cout<<sum<<"=1.0"<<endl;
		}
	}
}

//统计GP
void calculateGP(double** _gp, char** _mk, int sampleSize, int sampleIndex1, int sampleIndex2, double p[9][3], int ifPrint=VERBOSE_MODE) {
	double* gp = (double *) _gp;
	char* mk = (char*)_mk;
	/*
	 * 文件读出后比较BC1频率分布而得的频率数组gp
	 * Ethan：根据基因型统计对应的实际数据概率
	 */



	for (int i = 0; i < sampleSize; i++) {
		int lineNumber = 0;

		char cSample1 = *(mk+sampleIndex1*sampleSize+i);
		char cSample2 = *(mk+sampleIndex2*sampleSize+i);

		if (cSample1 == 'A') {
			if (cSample2 == 'A')
				lineNumber = 0;
			else if (cSample2 == 'H')
				lineNumber = 1;
			else if (cSample2 == 'B')
				lineNumber = 2;
		} else if (cSample1 == 'H') {
			if (cSample2 == 'A')
				lineNumber = 3;
			else if (cSample2 == 'H')
				lineNumber = 4;
			else if (cSample2 == 'B')
				lineNumber = 5;
		} else if (cSample1 == 'B') {
			if (cSample2 == 'A')
				lineNumber = 6;
			else if (cSample2 == 'H')
				lineNumber = 7;
			else if (cSample2 == 'B')
				lineNumber = 8;
		}

		*(gp+i*sampleSize+0) = p[lineNumber][0];
		*(gp+i*sampleSize+1) = p[lineNumber][1];
		*(gp+i*sampleSize+2) = p[lineNumber][2];
		//cout<<i <<"\t" <<cSample1 <<"\t" <<cSample2 <<"\t" <<lineNumber <<endl;
	}


	if (ifPrint) {
		cout<<"=======GPData======"<<endl;
		for (int c = 0; c < sampleSize; c++) {
			cout<<c <<"\t";
			for (int l = 0; l < 3; l++)
				cout<<*(gp+c*sampleSize+l) <<"\t";
			cout<<endl;
		}
	}
}

//计算LOD
double calculateLOD(double** _gp, int sampleSize, double* y, double s0, double s1, double u0, double u1, double u2, double u3, int ifPrint=VERBOSE_MODE) {
	double* gp = (double*) _gp;
	double sum1=0, sum2=0;
	//cout <<"u0=" <<u0 <<endl <<"s0=" <<s0 <<endl;
	for (int i = 0; i < sampleSize; i++) {
		sum1 += log10( f(y[i], u0, s0) );
		sum2 += log10( *(gp+i*sampleSize+0) * f(y[i], u1, s1)
				+ *(gp+i*sampleSize+1) * f(y[i], u2, s1)
				+ *(gp+i*sampleSize+2) * f(y[i], u3, s1)
			);
	}

	// LOD=log(L0/L)=logL0-LogL
	// 连乘Log改为加法，除法改减法
	double lod;
	lod = sum2 - sum1;

	if (ifPrint) {
		cout<<"sum1="<<sum1<<"\tsum2="<<sum2<<"\tLOD="<<lod<<endl;
	}
	return lod;
}

void EMCalculate(double* y, double** _gp, const int sampleSize, double u0, double s0, double* u1, double* u2, double* u3, double* s1, int ifPrint=VERBOSE_MODE) {
	double *gp = (double*)_gp;

	// u1, u2, u3, sigma2 见113
	double u10 = u0;
	double u20 = u0;
	double u30 = u0;
	double s10 = s0;//*s1+1;

	// 收敛精度要求
	int step=1;
	if (ifPrint) {
		cout<<"=======EMSteps======"<<endl;
		cout<<"[step]\tu1\tu2\tu3\ts1\n";
	}
	while ( fabs(u10 - *u1) > ACCURACY_REQ ||
			fabs(u20 - *u2) > ACCURACY_REQ ||
			fabs(u30 - *u3) > ACCURACY_REQ ||
			fabs(s10 - *s1) > ACCURACY_REQ) {
		double py1 = 0, py10 = 0;
		double py2 = 0, py20 = 0;
		double py3 = 0, py30 = 0;
		double py4 = 0, py5 = 0, py6 = 0;

		int i;
		for (i = 0; i < sampleSize; i++) {
			double gpi1f = *(gp+i*sampleSize+0) * f(y[i], u10, s10);
			double gpi2f = *(gp+i*sampleSize+1) * f(y[i], u20, s10);
			double gpi3f = *(gp+i*sampleSize+2) * f(y[i], u30, s10);
			double gpi = gpi1f + gpi2f + gpi3f;

			py1 += gpi1f / gpi;		//PI1的值
			py2 += gpi2f / gpi;	//PI2的值
			py3 += gpi3f / gpi;	//PI3的值

			py10 += gpi1f / gpi * y[i];	//PI1和PI1*YI的和
			py20 += gpi2f / gpi * y[i];	//PI2和PI2*YI的和
			py30 += gpi3f / gpi * y[i];	//PI3和PI3*YI的和

			py4 += gpi1f / gpi * (y[i] - u10) * (y[i] - u10);
			py5 += gpi2f / gpi * (y[i] - u20) * (y[i] - u20);
			py6 += gpi3f / gpi * (y[i] - u30) * (y[i] - u30);
		}

		*u1 = u10; u10 = py10 / py1;
		*u2 = u20; u20 = py20 / py2;
		*u3 = u30; u30 = py30 / py3;
		*s1 = s10; s10 = (py4 + py5 + py6) / sampleSize;

		if (ifPrint) {
			cout <<"[" <<step++ <<"]\t" <<*u1 <<"\t" <<*u2 <<"\t" <<*u3 <<"\t" <<*s1 <<endl;
		}
	}//end-while

	//输出终值
	*u1 = u10;
	*u2 = u20;
	*u3 = u30;
	*s1 = s10;
	if (ifPrint) {
		cout <<endl <<"Finals:" <<endl;
		cout <<"u1=" <<*u1 <<"\tu2=" <<*u2 <<"\tu3=" <<*u3 <<"\ts1=" <<*s1 <<endl;
	}
}

//TODO d2r是做什么的？
double d2r(double x) {
	return 0.5 * (1 - exp(-2 * x * 0.01) );
}

// 区间QTL分析
double intervalQTL(int currentTrait, double u0, double s0, double length, double startPoint, char** _mk, double* y, int sampleSize, int ifLast=0, int ifPrint=VERBOSE_MODE) {
	double r = d2r( length );
	double r1 = d2r ( startPoint );
	double r2 = d2r ( length-startPoint );

	if(ifPrint) {
		cout <<"r=" <<r <<"\tr1=" <<r1 <<"\tr2=" <<r2 <<"\td1=" <<startPoint <<"\td2=" <<(length-startPoint) <<"\td=" <<length <<endl;
	}

	//无干扰时F2群体 QTL 基因频率分布数组函数p 书107页//
	double p[9][3];
	calculatePij(p, r, r1, r2, ifPrint);

	//统计各类型分布比率
	double **gp = new double*[sampleSize];
	for (int i=0; i<sampleSize; i++) {
		gp[i] = new double[3];
	}

	calculateGP(gp, _mk, sampleSize, currentTrait, currentTrait+1, p, ifPrint);

	double u1=0, u2=0, u3=0, s1=0;
	// EM算法
	EMCalculate(y, gp, sampleSize, u0, s0, &u1, &u2, &u3, &s1, ifPrint);

	return calculateLOD(gp, sampleSize, y, s0, s1, u0, u1, u2, u3, ifPrint);
}

//打印帮助
void printHelp(char* fileName) {
	cout<<"Usage:\t" <<fileName <<" [ConfigFileName]" <<endl;
}

int main(int args, char* argv[]) {
	int iSampleSize;
	int iTraitNumber;
	double step;

	string sExpressDataFile;
	string sGeneDataFile;
	string sTraitIntervalFile;

	bool ifPrintInitDataReport;
	bool ifPrintCalReport;
	bool ifPrintFinalReport;

	if (args==1) {
		cout<<"Please input the Number of Samples:";
		cin>>iSampleSize;

		cout<<iSampleSize<<endl<<"Please input the Number of Traits: ";
		cin>>iTraitNumber;

		cout<<iTraitNumber<<endl<<"Please input the length of step: ";
		cin>>step;

		cout<<step<<endl<<"Please locate the file of Express Data: ";
		cin>>sExpressDataFile;

		cout<<sExpressDataFile<<endl<<"Please locate the file of Gene Data: ";
		cin>>sGeneDataFile;

		cout<<sGeneDataFile<<endl<<"Please locate the file of Trait Interval Data: ";
		cin>>sTraitIntervalFile;

		cout<<sTraitIntervalFile<<endl<<"Do you want detailed Initialization Report?(0/1) ";
		cin>>ifPrintInitDataReport;

		cout<<(ifPrintInitDataReport?"True":"False")<<endl<<"Do you want detailed Calculation Report?(0/1) ";
		cin>>ifPrintCalReport;

		cout<<(ifPrintCalReport?"True":"False")<<endl<<"Do you want detailed Final Report?(0/1) ";
		cin>>ifPrintFinalReport;

		cout<<(ifPrintFinalReport?"True":"False")<<endl;
	} else if (args==2){
		cout <<"Input File:" <<argv[1] <<endl;
		fstream fin;
		fin.open(argv[1],ios::in);

		fin >>iSampleSize >>iTraitNumber >>step;
		fin >>sExpressDataFile;
		fin >>sGeneDataFile;
		fin >>sTraitIntervalFile;
		fin >>ifPrintInitDataReport >>ifPrintCalReport >>ifPrintFinalReport;

		fin.close();
	} else {
		printHelp(argv[0]);
	}

	cout <<"Input Summary:\t" <<iSampleSize <<" samples of " <<iTraitNumber <<" traits in following files:" <<endl;
	cout <<"Express Data in:        " <<sExpressDataFile <<endl;
	cout <<"Gene Data in:           " <<sGeneDataFile <<endl;
	cout <<"Trait Interval Data in: " <<sTraitIntervalFile <<endl;
	cout <<"Log Detail: " <<endl
			<<"\tInitialization Report - " <<(ifPrintInitDataReport?"Yes":"No") <<endl
			<<"\tCalculation Report - " <<(ifPrintCalReport?"Yes":"No") <<endl
			<<"\tFinal Report - " <<(ifPrintFinalReport?"Yes":"No") <<endl;
	cout <<"Calculate step:\t" <<step <<endl;

	// 表型数据
	double* y = new double[iSampleSize];
	// 读取表型数据文件，并计算表型数组的均值、方差
	double u0=0, s0=0;
	initExpressData(sExpressDataFile, iSampleSize, y, &u0, &s0, ifPrintInitDataReport);

	// 基因型数据
	//char mk[TRAIT_NUMBER][SAMPLE_SIZE];
	char** mk = new char*[iTraitNumber];
	for (int i=0; i<iTraitNumber; i++) {
		mk[i] = new char[iSampleSize];
	}

	//读取基因型数据文件
	initGeneData(sGeneDataFile, iSampleSize, iTraitNumber, mk, ifPrintInitDataReport);

	//读取位点距离数据，来源为见书115页
	double* traitInterval= new double[iTraitNumber];
	initIntervalData(sTraitIntervalFile, iTraitNumber-1, traitInterval, ifPrintInitDataReport);

	//逐个位点计算LOD值
	for (int currentTrait = 0; currentTrait<iTraitNumber-1; currentTrait++) {
		if (ifPrintFinalReport) {
			cout<<"---------"<<(currentTrait+1)<<"---("<<traitInterval[currentTrait]<<")--------"<<endl;
		}
		for (double startPoint=0.0; startPoint<traitInterval[currentTrait]; startPoint += step) {
			double LOD = intervalQTL(currentTrait, u0, s0, traitInterval[currentTrait], startPoint, mk, y, iSampleSize, currentTrait==iTraitNumber-1, ifPrintCalReport);
			if (ifPrintFinalReport) {
				cout<<"["<<startPoint<<"~"<<(startPoint+step)<<"]:["<<traitInterval[currentTrait]<<"] LOD="<<LOD<<endl;
			}
		}
	}
	return 1;
}
