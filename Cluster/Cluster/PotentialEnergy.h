#pragma once
#include "Atom.h"
#include "Clusters.h"

enum PE_TYPE
{
	PE_LJ,
	PE_Gupta,
	PE_FS,
	PE_SC
};
struct PE_AtomParamter{ };

class PotentialEnergy
{
public:
	static PotentialEnergy* initWithType(PE_TYPE type);

	virtual ~PotentialEnergy();
	
	virtual double EnergyValue(Clusters& cluster);

	virtual double ForceValue(Clusters& cluster);

//	virtual double atomEnergyAtIndex(int index);

private:
	virtual PE_AtomParamter& ReturnAtomParameter(ATOM_TYPE atom);

	virtual PE_AtomParamter& ReturnAtomParameter(ATOM_TYPE atom1,ATOM_TYPE atom2);


protected:
	vector<double> _atomEnergy;
	vector<double> _atomForceAlongX;
	vector<double> _atomForceAlongY;
	vector<double> _atomForceAlongZ;
};

string PE_EnergyTypeName( PE_TYPE type );