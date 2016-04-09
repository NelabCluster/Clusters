#pragma once
#include "clusters.h"

class DEBaseClusters :
	public Clusters
{
public:
	DEBaseClusters(int N,vector<ATOM_TYPE> types,vector<int> numbers);
	~DEBaseClusters(void);

	friend DEBaseClusters operator+(DEBaseClusters& left,DEBaseClusters& right);

	friend DEBaseClusters operator-(DEBaseClusters& left,DEBaseClusters& right);

	friend DEBaseClusters operator*(double value,DEBaseClusters& right);

	friend DEBaseClusters operator*(DEBaseClusters& left,double value);
};

