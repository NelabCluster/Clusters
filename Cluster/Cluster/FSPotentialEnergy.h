#pragma once
#include "potentialenergy.h"

struct FS_AtomParamter : PE_AtomParamter
{
	double d;
	double A;
	double beta;
	double c;
	double c0;
	double c1;
	double c2;

	FS_AtomParamter(double d,double A,double beta,double c,double c0,double c1,double c2)
		:d(d),A(A),beta(beta),c(c),c0(c0),c1(c1),c2(c2){ }
};

static FS_AtomParamter FS_Fe(3.569745,1.828905,1.8,3.40,1.2371147,-0.3592185,-0.0385607);

class FSPotentialEnergy :
	public PotentialEnergy
{
public:
	FSPotentialEnergy(void);
	~FSPotentialEnergy(void);
	
	double EnergyValue(Clusters& cluster);
	double ForceValue(Clusters& cluster);

private:
	PE_AtomParamter& ReturnAtomParameter(ATOM_TYPE atom);

	PE_AtomParamter& ReturnAtomParameter(ATOM_TYPE atom1,ATOM_TYPE atom2);
};

