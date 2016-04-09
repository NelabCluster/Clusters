#include "StdAfx.h"
#include "PotentialEnergy.h"
#include "LJPotentialEnergy.h"
#include "GuptaPotentialEnergy.h"
#include "FSPotentialEnergy.h"

PotentialEnergy* PotentialEnergy::initWithType(PE_TYPE type)
{
	switch(type)
	{
	case PE_LJ:
		return new LJPotentialEnergy();
	case PE_Gupta:
		return new GuptaPotentialEnergy();
	case PE_FS:
		return new FSPotentialEnergy();
	}
	return new PotentialEnergy();
}
PotentialEnergy::~PotentialEnergy(){ }

double PotentialEnergy::EnergyValue(Clusters& cluster){ return 0; }

double PotentialEnergy::ForceValue(Clusters& cluster){ return 0; }

PE_AtomParamter& PotentialEnergy::ReturnAtomParameter(ATOM_TYPE atom){ return PE_AtomParamter(); }

PE_AtomParamter& PotentialEnergy::ReturnAtomParameter(ATOM_TYPE atom1,ATOM_TYPE atom2){ return PE_AtomParamter(); }

string PE_EnergyTypeName( PE_TYPE type )
{
	switch( type )
	{
		case PE_LJ: return "LJ";
		case PE_Gupta: return "Gupta";
		case PE_FS: return "FS";
		case PE_SC: return "SC";
	}
}