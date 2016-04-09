#pragma once
#include "Common.h"

class InfoModel
{
public:
	InfoModel(void);
	~InfoModel(void);

	double energy;

	int N;

	Alloy alloy;

	AlloyNum alloyNum;

	static void ModelWithCluster( InfoModel &model, Clusters &cluster);
	static bool ModelFromFile( InfoModel& model, string &file );

	void SaveToFile( string &fileName );
};

