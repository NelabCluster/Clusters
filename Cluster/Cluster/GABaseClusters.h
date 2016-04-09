#pragma once
#include "clusters.h"
class GABaseClusters :
	public Clusters
{
public:
	GABaseClusters(int N,vector<ATOM_TYPE> types,vector<int> numbers);
	~GABaseClusters(void);

};

