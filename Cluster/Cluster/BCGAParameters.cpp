#include "StdAfx.h"
#include "BCGAParameters.h"


BCGAParameters::BCGAParameters(int N,PE_TYPE energy, ATOM_TYPE type,string root)
	:AlgorithmParameters(N,energy,type,root) { this->defaultPara(); }

BCGAParameters::BCGAParameters(int N, PE_TYPE energyType,ATOM_TYPE atom1,ATOM_TYPE atom2,int numberOfAtom1,string root)
	:AlgorithmParameters(N,energyType,atom1,atom2,numberOfAtom1,root) { this->defaultPara(); }


BCGAParameters::~BCGAParameters(void) { }

void BCGAParameters::defaultPara()
{
	popSize = 20;
	gen = 50;
	crossoverRate = 0.8;
	muteRate = 0.3;
	init = BCGA_Init_Swap;
	crossover = BCGA_Crossover_PlaneCutSplice;
	mutate = BCGA_Mutate_Displace;
	mutateNum = 1;

	if ( atomTypes.size() >=2 )
		mutate |= BCGA_Mutate_Swap;
}
