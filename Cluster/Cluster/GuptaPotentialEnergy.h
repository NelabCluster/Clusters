#pragma once
#include "potentialenergy.h"

struct Gupta_AtomParamter : PE_AtomParamter
{
	double r0;
	double A;
	double Xi;
	double P;
	double q;
	Gupta_AtomParamter(double r0,double A,double Xi,double P,double q)
		:r0(r0),A(A),Xi(Xi),P(P),q(q){ }
};

static Gupta_AtomParamter GuptaFe(2.553,0.13315,1.6179,10.5,2.6),
						  GuptaPt(2.7746,0.2975,2.695,10.612,4.004),
				          GuptaFePt(2.6638,0.19903,2.0881,10.556,3.302),
					      //IntJQuantumChem.25(2003)41
						  GuptaCu(2.556,0.0855,1.2240,10.960,2.2780),
						  GuptaAu(2.884,0.2061,1.7900,10.229,4.0360),
						  GuptaCuAu(2.556,0.1539,1.5605,11.050,3.0475),
						  //J. Chem. Phys. 122, 244707 (2005)
						  GuptaCo(2.497,0.0950,1.4880,11.604,2.286),
						  //PhysRevB.68(2003)195418 论文里将r0归一化
						  GuptaZn(1,0.1477,0.8900,9.689,4.602),		//r0 = 2.66
						  GuptaCd(1.0,0.1420,0.8117,10.612,5.206);	//r0 = 2.98

class GuptaPotentialEnergy :
	public PotentialEnergy
{
public:
	GuptaPotentialEnergy(void);
	~GuptaPotentialEnergy(void);


	double EnergyValue(Clusters& cluster);

	double ForceValue(Clusters& cluster);

private:
	PE_AtomParamter& ReturnAtomParameter(ATOM_TYPE atom);

	PE_AtomParamter& ReturnAtomParameter(ATOM_TYPE atom1,ATOM_TYPE atom2);
};

