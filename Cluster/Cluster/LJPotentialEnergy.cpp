#include "StdAfx.h"
#include "LJPotentialEnergy.h"


LJPotentialEnergy::LJPotentialEnergy(void)
{
}


LJPotentialEnergy::~LJPotentialEnergy(void)
{
}

double LJPotentialEnergy::EnergyValue(Clusters& cluster)
{
	int N = cluster.GetAtomsNumber();
	double *dis = cluster.GetDistancePointer();
	double E = 0;
	vector<double> EOfAtom( N, 0 );
	
	for (int i = 0; i < N - 1; i ++)
	{
		for (int j = i + 1; j < N; j ++)
		{
			double r = dis[ i * N + j ];
			double e = pow(1/r,12) - pow(1/r,6);
			EOfAtom[ i ] += 2 * e;
			EOfAtom[ j ] += 2 * e; 
			E += 4 * e;
		}
	}

	cluster.SetEnergyVectorOfAtoms( EOfAtom );
	cluster.SetEnergy( E );

	return E;
}

double LJPotentialEnergy::ForceValue( Clusters& cluster )
{
	int N = cluster.GetAtomsNumber();
	double *dis = cluster.GetDistancePointer();
	vector<double> FAlongX(N,0), FAlongY(N,0), FAlongZ(N,0);

	for (int i = 0; i < N - 1; i ++)
	{
		for (int j = i + 1; j < N; j ++)
		{
			double r = dis[ i * N + j ];
			double FK = 4 * ( -12 * pow(1/r,13) + 6 * pow(1/r,7) );

			AtomPos posOfAtomAtI = cluster.GetAtomAtIndex(i).GetPos();
			AtomPos posOfAtomAtJ = cluster.GetAtomAtIndex(j).GetPos();
			FAlongX[i] = FAlongX[i]+FK*(posOfAtomAtI.x-posOfAtomAtJ.x)/r;
			FAlongX[j] = FAlongX[j]-FK*(posOfAtomAtI.x-posOfAtomAtJ.x)/r;
			FAlongY[i] = FAlongY[i]+FK*(posOfAtomAtI.y-posOfAtomAtJ.y)/r;
			FAlongY[j] = FAlongY[j]-FK*(posOfAtomAtI.y-posOfAtomAtJ.y)/r;
			FAlongZ[i] = FAlongZ[i]+FK*(posOfAtomAtI.z-posOfAtomAtJ.z)/r;
			FAlongZ[j] = FAlongZ[j]-FK*(posOfAtomAtI.z-posOfAtomAtJ.z)/r;
		}
	}

	cluster.SetForceXYZVectorOfAtoms(FAlongX,FAlongY,FAlongZ);

	double forceXOfCluster = 0, forceYOfCluster = 0, forceZOfCluster = 0;
	for (int i = 0; i < N; i++)
	{
		forceXOfCluster += FAlongX[i] * FAlongX[i];
		forceYOfCluster += FAlongY[i] * FAlongY[i];
		forceZOfCluster += FAlongZ[i] * FAlongZ[i];
	}

	double F = sqrt(forceXOfCluster  + forceYOfCluster  + forceZOfCluster );
	cluster.SetForceOfCluster(F);

	return F;
}