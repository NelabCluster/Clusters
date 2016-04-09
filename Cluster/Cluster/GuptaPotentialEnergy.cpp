#include "StdAfx.h"
#include "GuptaPotentialEnergy.h"


GuptaPotentialEnergy::GuptaPotentialEnergy(void)
{
}


GuptaPotentialEnergy::~GuptaPotentialEnergy(void)
{

}


double GuptaPotentialEnergy::EnergyValue(Clusters& cluster)
{
	int N = cluster.GetAtomsNumber();
	_atomEnergy.resize(N);
	double *dis = cluster.GetDistancePointer();
	double E = 0;
	vector<double> VEN(N,0);
	vector<double> PEN(N,0);
	
	Alloy alloy = cluster.GetAlloy();

	for (int i = 0; i < N - 1; i ++)
	{
		for (int j = i + 1; j < N; j ++)
		{
			double r = dis[ i * N + j ];
			int note1 = cluster.GetAtomAtIndex(i).GetNote();
			int note2 = cluster.GetAtomAtIndex(j).GetNote();
			Gupta_AtomParamter parameter = (Gupta_AtomParamter&)ReturnAtomParameter(alloy[note1],alloy[note2]);
			double FMJN = r / parameter.r0 - 1;
			double FMJV = parameter.A * exp( -parameter.P * FMJN );
			double FMJP  = parameter.Xi * parameter.Xi * exp( -2 * parameter.q *FMJN );

			VEN[i] += FMJV;
			VEN[j] += FMJV;
			PEN[i] += FMJP;
			PEN[j] += FMJP;
		}
	}
	for (int i = 0; i < N; i++)
	{		
		_atomEnergy[i] = VEN[i] - sqrt(PEN[i]);
		E += _atomEnergy[i];
	}

	cluster.SetEnergyVectorOfAtoms(_atomEnergy);
	cluster.SetEnergy(E);

	return E;
}

double GuptaPotentialEnergy::ForceValue(Clusters& cluster)
{
	int N = cluster.GetAtomsNumber();
	double *dis = cluster.GetDistancePointer();
	vector<double> PEN(N,0);
	Alloy alloy = cluster.GetAlloy();

	for (int i = 0; i < N - 1; i ++)
	{
		for (int j = i + 1; j < N; j ++)
		{
			double r = dis[ i * N + j ];
			int note1 = cluster.GetAtomAtIndex(i).GetNote();
			int note2 = cluster.GetAtomAtIndex(j).GetNote();
			Gupta_AtomParamter parameter = (Gupta_AtomParamter&)ReturnAtomParameter(alloy[note1],alloy[note2]);
			double FMJN = r / parameter.r0 - 1;
			double FMJP  = parameter.Xi * parameter.Xi * exp( -2 * parameter.q *FMJN );

			PEN[i] += FMJP;
			PEN[j] += FMJP;
		}
	}

	for (int i = 0; i < N; i++)
		PEN[i] = 1/sqrt(PEN[i])/2;

	_atomForceAlongX.resize(N);
	_atomForceAlongY.resize(N);
	_atomForceAlongZ.resize(N);
	_atomForceAlongX.assign(N,0);
	_atomForceAlongY.assign(N,0);
	_atomForceAlongZ.assign(N,0);

	for (int i = 0; i < N-1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			double r = dis[ i * N + j ];
			int note1 = cluster.GetAtomAtIndex(i).GetNote();
			int note2 = cluster.GetAtomAtIndex(j).GetNote();
			Gupta_AtomParamter parameter = (Gupta_AtomParamter&)ReturnAtomParameter(alloy[note1],alloy[note2]);
			double FMJN = r / parameter.r0 - 1;

			double dFMJV = - parameter.P / parameter.r0 * parameter.A * exp( -parameter.P * FMJN ) * 2;
			double dFMJP = -2 * parameter.q / parameter.r0 * parameter.Xi * parameter.Xi * exp( -2 * parameter.q *FMJN );
			double FK = dFMJV - (PEN[i] + PEN[j]) * dFMJP;

			Atom atomOfI = cluster.GetAtomAtIndex(i);
			AtomPos posOfAtomAtI = atomOfI.GetPos();
			Atom atomOfJ = cluster.GetAtomAtIndex(j);
			AtomPos posOfAtomAtJ = atomOfJ.GetPos();
			_atomForceAlongX[i] = _atomForceAlongX[i]+FK*(posOfAtomAtI.x-posOfAtomAtJ.x)/r;
			_atomForceAlongX[j] = _atomForceAlongX[j]-FK*(posOfAtomAtI.x-posOfAtomAtJ.x)/r;
			_atomForceAlongY[i] = _atomForceAlongY[i]+FK*(posOfAtomAtI.y-posOfAtomAtJ.y)/r;
			_atomForceAlongY[j] = _atomForceAlongY[j]-FK*(posOfAtomAtI.y-posOfAtomAtJ.y)/r;
			_atomForceAlongZ[i] = _atomForceAlongZ[i]+FK*(posOfAtomAtI.z-posOfAtomAtJ.z)/r;
			_atomForceAlongZ[j] = _atomForceAlongZ[j]-FK*(posOfAtomAtI.z-posOfAtomAtJ.z)/r;
		}
	}

	cluster.SetForceXYZVectorOfAtoms(_atomForceAlongX,_atomForceAlongY,_atomForceAlongZ);

	double forceXOfCluster = 0, forceYOfCluster = 0, forceZOfCluster = 0;
	for (int i = 0; i < N; i++)
	{
		forceXOfCluster += _atomForceAlongX[i] * _atomForceAlongX[i];
		forceYOfCluster += _atomForceAlongY[i] * _atomForceAlongY[i];
		forceZOfCluster += _atomForceAlongZ[i] * _atomForceAlongZ[i];
	}

	double F = sqrt(forceXOfCluster  + forceYOfCluster  + forceZOfCluster );
	cluster.SetForceOfCluster(F);

	return F;
}

PE_AtomParamter& GuptaPotentialEnergy::ReturnAtomParameter(ATOM_TYPE atom){ 
	
	switch (atom)
	{
	case ATOM_Fe:
		return GuptaFe;
		break;
	case ATOM_Pt:
		return GuptaPt;
	case ATOM_Cu:
		return GuptaCu;
	case ATOM_Au:
		return GuptaAu;
	case ATOM_Co:
		return GuptaCo;
	case ATOM_Zn:
		return GuptaZn;
	case ATOM_Cd:
		return GuptaCd;
		break;
	}

	return Gupta_AtomParamter(0,0,0,0,0);
}

PE_AtomParamter& GuptaPotentialEnergy::ReturnAtomParameter(ATOM_TYPE atom1,ATOM_TYPE atom2)
{
	if ( atom1 == atom2)
	{
		return ReturnAtomParameter(atom1);
	}

	if ( (atom1 == ATOM_Fe && atom2 == ATOM_Pt) || (atom1 == ATOM_Pt && atom2 == ATOM_Fe))
	{
		return GuptaFePt;
	}
	if ( (atom1 == ATOM_Cu && atom2 == ATOM_Au) || (atom1 == ATOM_Au && atom2 == ATOM_Cu))
	{
		return GuptaCuAu;
	}

	return Gupta_AtomParamter(0,0,0,0,0);
}