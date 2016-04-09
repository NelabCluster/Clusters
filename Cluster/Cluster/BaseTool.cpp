#include "StdAfx.h"
#include "BaseTool.h"


BaseTool::BaseTool(void)
{
}


BaseTool::~BaseTool(void)
{
}

vector<int> BaseTool::RandPerm(int N, int K)
{
	vector<int> allP(N);
	vector<int> p(K);

	for(int i = 0; i < N; i++)
		allP[i] = i;
	for(int i=0;i < K;i++)
	{
		int point = RANDIUINT(N-i);
		int temp = allP[i];
		allP[i] = allP[point+i];
		allP[point+i] = temp;
		p[i] = allP[i];
	}
	return p;
}

string BaseTool::IntToString( int m )
{
	stringstream ss;
	ss<<m;
	return ss.str();
}

string BaseTool::DoubleToString( double m, int decimal )
{
	stringstream ss;
	ss<<setprecision(decimal)<<m;
	return ss.str();
}