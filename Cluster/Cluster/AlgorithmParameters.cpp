#include "StdAfx.h"
#include "AlgorithmParameters.h"

AlgorithmParameters::AlgorithmParameters( int n, PE_TYPE pe, ATOM_TYPE atomType, string root )
{
	N = n;
	energyType = pe;
	atomTypes.push_back(atomType);
	atomNumber.push_back(N);
	rootFolderName = root;
}

AlgorithmParameters::AlgorithmParameters(int n,PE_TYPE pe,ATOM_TYPE atom1,ATOM_TYPE atom2,int atom1N, string root)
{
	N = n;

	energyType = pe;

	if ( NULL != atom1 ) atomTypes.push_back(atom1);

	if ( NULL != atom2 ) atomTypes.push_back(atom2);

	atomNumber.push_back(atom1N);

	atomNumber.push_back(N-atom1N);

	rootFolderName = root;
}

AlgorithmParameters::AlgorithmParameters(int n, PE_TYPE pe, Alloy &types, AlloyNum &numbers, string root)
{
	N = n;
	energyType = pe;
	atomTypes = types;
	atomNumber = numbers;
	rootFolderName = root;
}


AlgorithmParameters::~AlgorithmParameters(void)
{

}
