#pragma once
#include "algorithmparameters.h"


enum BCGA_Init
{
	BCGA_Init_Rand,
	BCGA_Init_Grow,
	BCGA_Init_Swap
};

enum BCGA_Crossover
{
	BCGA_Crossover_Value,
	BCGA_Crossover_PlaneCutSplice
};

enum BCGA_Mutate
{
	BCGA_Mutate_Displace = 0x01,
	BCGA_Mutate_Twitst = 0x02,
	BCGA_Mutate_Replace = 0x04,
	BCGA_Mutate_Swap = 0x08
};

class BCGAParameters :
	public AlgorithmParameters
{
public:
	BCGAParameters( int N, PE_TYPE energy, ATOM_TYPE type, string root);
	BCGAParameters( int N, PE_TYPE energyType, ATOM_TYPE atom1, ATOM_TYPE atom2, int numberOfAtom1, string root);

	~BCGAParameters(void);

	int popSize;

	double crossoverRate;

	double muteRate;

	int gen;
	
	BCGA_Init init;

	BCGA_Crossover crossover;

	int mutate;

	int mutateNum;

private:
	void defaultPara();
};

