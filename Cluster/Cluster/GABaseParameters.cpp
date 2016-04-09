#include "StdAfx.h"
#include "GABaseParameters.h"


GABaseParameters::GABaseParameters(int N, PE_TYPE energyType,ATOM_TYPE atom1,ATOM_TYPE atom2,int numberOfAtom1,string root)
	:AlgorithmParameters(N,energyType,atom1,atom2,numberOfAtom1,root)
{
	
}


GABaseParameters::~GABaseParameters(void)
{
}
