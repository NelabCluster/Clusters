#pragma once
#include "Common.h"
#include "FileTool.h"

class ClusterPool
{

public:

	static ClusterPool* Share() 
	{
		static ClusterPool* instance;
		if (instance == nullptr)
			instance = new ClusterPool();
		return instance;
	}

	void SaveCluster( Clusters &cluster, PE_TYPE energyType );

	double GetEnergy( PE_TYPE energyType, Alloy &alloy, AlloyNum &num, int N );
	void SetInfo( Clusters &cluster, PE_TYPE energyType );

	bool GetCluster( Clusters &cluster, PE_TYPE energyType );
	void SetCluster( Clusters &cluster, PE_TYPE energyType );

	void SetDiamomd( Clusters &cluster, PE_TYPE energyType  );

private:


};
