#pragma once
#include "IAlgorithm.h"

class AlgorithmManage
{
public:
	AlgorithmManage( IAlgorithm* a );
	AlgorithmManage(void);
	~AlgorithmManage(void);

	void start();
	void start( int repeat );

private:
	IAlgorithm* _algorithm;

	void _Start();
};

