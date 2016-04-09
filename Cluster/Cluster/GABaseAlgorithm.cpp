#include "StdAfx.h"
#include "GABaseAlgorithm.h"




GABaseAlgorithm::GABaseAlgorithm(GABaseParameters& para)
	:_parameters(para)
{

}

GABaseAlgorithm::~GABaseAlgorithm(void)
{
}

void GABaseAlgorithm::Initialization()
{
	double cubic = pow(_parameters.N*1.0,1/3);
	for (int i =0 ;i < _parameters._numberOfClusters; i++)
	{
		GABaseClusters newCluster(_parameters.N,_parameters.atomTypes,_parameters.atomNumber);
		//newCluster.RandInCubic(cubic);
		//newCluster.GetEnergyWithLocal(_parameters._energyType);
		_pop.push_back(newCluster);
	}
}

void GABaseAlgorithm::Process()
{

}

bool GABaseAlgorithm::EndCondition()
{
	return true;
}