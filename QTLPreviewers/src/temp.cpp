#include "IMQTLPreviewer.h"

template<class Iterator>
void Myrandom_shuffle(Iterator beg,Iterator end) {
	Iterator it;
	srand(unsigned(time(0)));
	int diff;
    Iterator swp;

	for(it=beg+1;it!=end;it++) {
		diff=rand()%(it-beg+1);
		swp=beg+diff;
		swap(*it,*swp);
	}
}

typedef int Type;

template<class Iterator>
inline static void Output(Iterator beg,Iterator end) {
	copy(beg,end,ostream_iterator<Type>(cout," "));
}

void testR() {
	int n=300;
	float t;
	vector<Type>v;
	int i;
	cout<<"Original:"<<endl;
	for(i=0;i<n;i++) {
		v.push_back(i+1);
	}
	Output(v.begin(),v.end());

	cout<<endl;

	getTimeStamp();
	int times = 1000;
	for (int i=0; i<times; i++) {
		Myrandom_shuffle(v.begin(),v.end());
	}
	t = getTimeStamp();

	cout <<"Total time: " <<t <<" ms" <<endl
			<<"Average time:" << t/times <<" ms" <<endl;

	cout<<"Reordered:"<<endl;
	Output(v.begin(),v.end());
	cout<<endl;
}

template<class Iterator>
double getAverageNew(Iterator begin, Iterator end) {
	Iterator it;
	double sum = 0;
	for(it=begin;it!=end;it++) {
		sum += *it;
	}
	return sum/double(end-begin);
}

int tempMain() {
	vector<EXPData> expData;
	fstream fin;
	fin.open("data/express-66.txt", ios::in);
	while (!fin.eof()) {
		EXPData e;
		fin >> e;
		expData.push_back(e);
	}

	EXPData u0 = getAverageNew( expData.begin(), expData.end() );
	EXPData s0 = getVariance(expData.begin(), expData.end(), u0);

	/*
	for (int i = 0; i < 189; i++) {
		EXPData input;
		fin >> data[i];
		//data.push_back(s);
	}
	*/
	fin.close();

	cout << "=======ExpressData=======" << endl;
	for (int i = 0; i < 189; i++) {
		cout << "[" << i << "]\t" << expData[i] << endl;
	}
	cout << "u0=" << u0 << endl << "s0=" << s0 << endl;

	return 0;
}
