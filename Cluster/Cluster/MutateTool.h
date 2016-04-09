#pragma once
#include "Clusters.h"

class MutateTool
{
public:

	static void AtomDisplacement(Clusters& cluster,int number);
	
	static void Twisting(Clusters& cluster);

	static void Replacement(Clusters& cluster);

	static void AtomRandSwap(Clusters& cluster, int times);
};

