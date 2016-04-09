#pragma once
#include "Common.h"


class InitTool
{
public:

	static void RandInCubic( Clusters& cluster, double length );
	static void RandInSphere( Clusters& cluster, double radius );
	static void GrowFromSmallerCluster( Clusters& cluster, Clusters& smallClusters );
	
	static void AtomTypeSwap( Clusters &cluster, Clusters &neighCluster, PE_TYPE type  ); 
	static void RandType( Clusters &cluster );
};

