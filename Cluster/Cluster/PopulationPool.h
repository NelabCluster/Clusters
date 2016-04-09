#pragma once
#include "Common.h"
class PopulationPool
{
public:
	PopulationPool(void);
	~PopulationPool(void);

private:
	vector<Clusters> pop;
};

