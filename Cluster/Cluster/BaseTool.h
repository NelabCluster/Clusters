#pragma once
class BaseTool
{
public:
	BaseTool(void);
	~BaseTool(void);

	static vector<int> RandPerm(int N,int K);
	static string IntToString( int m );
	static string DoubleToString( double m, int decimal );
};

