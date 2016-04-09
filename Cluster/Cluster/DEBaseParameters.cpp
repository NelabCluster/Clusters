#include "StdAfx.h"
#include "DEBaseParameters.h"


DEBaseParameters::DEBaseParameters(int N, PE_TYPE energyType,ATOM_TYPE atom1,ATOM_TYPE atom2,int numberOfAtom1,string root)
	:AlgorithmParameters(N,energyType,atom1,atom2,numberOfAtom1,root)
{
	_mute = DERAND1;
	_numberOfClusters = N + 10;
	_maxIteration = 200;
}


DEBaseParameters::~DEBaseParameters(void)
{
}
