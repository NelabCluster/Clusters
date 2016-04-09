#pragma once
#include "Common.h"
#include "PotentialEnergy.h"

enum Local_Type {
	LOCALTYPE_LBFGS
};

class LocalTool
{
public:
	LocalTool(void);
	~LocalTool(void);

	static double LocalMinimize(Clusters& beforeLocal, Clusters& afterLocal,PE_TYPE type);
	static double LocalMinimizeClusters(Clusters& cluster,PE_TYPE type);
};

