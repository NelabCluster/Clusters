#pragma once
#include "PotentialEnergy.h"
#include "Common.h"

class AlgorithmParameters
{
public:
	AlgorithmParameters( int n, PE_TYPE pe, ATOM_TYPE atomType, string root );
	AlgorithmParameters( int n, PE_TYPE pe, ATOM_TYPE atom1, ATOM_TYPE atom2, int atom1N, string root );
	AlgorithmParameters( int n, PE_TYPE pe, Alloy &types, AlloyNum &numbers, string root );
	~AlgorithmParameters(void);

	int N;

	PE_TYPE energyType;

	Alloy atomTypes;

	AlloyNum atomNumber;

	string rootFolderName;
};

