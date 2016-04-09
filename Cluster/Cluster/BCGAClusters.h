#pragma once
#include "Clusters.h"

class BCGAClusters:
	public Clusters
{
public:
	BCGAClusters(int N,vector<ATOM_TYPE> types,vector<int> numbers);
	~BCGAClusters(void);

	double fitness;
};

