#pragma once
#include "algorithmparameters.h"
class GABaseParameters :
	public AlgorithmParameters
{
public:
	GABaseParameters(int N, PE_TYPE energyType,
					  ATOM_TYPE atom1,ATOM_TYPE atom2,int numberOfAtom1,string root);
	~GABaseParameters(void);

	int _numberOfClusters;

	double _crossoverRate;

	double _muteRate;
};

