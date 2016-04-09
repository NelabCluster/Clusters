#include "StdAfx.h"
#include "FSPotentialEnergy.h"


FSPotentialEnergy::FSPotentialEnergy(void)
{
}


FSPotentialEnergy::~FSPotentialEnergy(void)
{
}

double FSPotentialEnergy::EnergyValue(Clusters& cluster)
{
	int i,j;
	double d,A,beta,c,c0,c1,c2;
	double r;
	double tempV,tempP;
	double *VEN,*PEN;
	double E = 0;
	double minR = 1.6,tempPmin;
	int N;
	double *R;
	//Fe
	d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;

	N = cluster.GetAtomsNumber();
	VEN = (double *)calloc(N,sizeof(double));
	PEN = (double *)calloc(N,sizeof(double));
	R = cluster.GetDistancePointer();

	tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
	for(i = 0; i < N-1; i++)
	{
		for(j = i+1; j < N; j++)
		{
			r = *(R + i*N + j);
			tempV = (r>=c)?0 : (r-c)*(r-c)*(c0+c1*r+c2*r*r);
			tempP = (r>=d)?0 : ((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);

			VEN[i] += tempV/2;
			VEN[j] += tempV/2;
			PEN[i] += tempP;
			PEN[j] += tempP;
		}
	}

	for(i=0 ;i < N; i++){
		PEN[i] = sqrt(PEN[i]);
		E += VEN[i] - A * PEN[i];
	}

	cluster.SetEnergy(E);

	free(VEN);
	free(PEN);
	return E;
}

double FSPotentialEnergy::ForceValue(Clusters& cluster)
{
	int i,j;	
	double d,A,beta,c,c0,c1,c2;
	double r;
	double tempdV,tempP,tempdP;
	double *PEN;
	double maxF = 0,FK;
	double minR = 1.6,tempPmin;
	double *FX,*FY,*FZ;
	int N;
	double *R;

	N = cluster.GetAtomsNumber();
	PEN = (double *)calloc(N,sizeof(double));
	R = cluster.GetDistancePointer();

	FX = (double *)calloc(N,sizeof(double)); memset(FX,0,N * sizeof(double));
	FY = (double *)calloc(N,sizeof(double)); memset(FY,0,N * sizeof(double));
	FZ = (double *)calloc(N,sizeof(double)); memset(FY,0,N * sizeof(double));

	//Fe
	d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;

	tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
	for(i=0;i<N-1; i++)
	{
		for(j=i+1;j<N;j++)
		{
			r = *(R + i*N + j);
			tempP = (r>=d)?0:((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);

			PEN[i] += tempP;
			PEN[j] += tempP;
		}
	}

	for(i = 0; i < N ;i ++)
		PEN[i] = (PEN[i] == 0)?0:1 / sqrt(PEN[i]) / 2;

	for(i = 0; i < N-1;i ++)
	{
		for(j=i+1;j<N;j++)
		{
			r = *(R + i*N + j);

			tempdV = (r>=c)?0 : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
			tempdP = (r>=d)? 0 : (2 * (r -d) + 3 * beta * (r - d) * (r - d) / d);

			tempdP = (PEN[i] + PEN[j]) * tempdP;

			FK = -tempdV / 2 + A * tempdP;

			AtomPos atomOfIPos = cluster.GetAtomAtIndex(i).GetPos();
			AtomPos atomOfJPos = cluster.GetAtomAtIndex(j).GetPos();
			FX[i]=FX[i]+FK*(atomOfIPos.x-atomOfJPos.x)/r;
			FX[j]=FX[j]-FK*(atomOfIPos.x-atomOfJPos.x)/r;
			FY[i]=FY[i]+FK*(atomOfIPos.y-atomOfJPos.y)/r;
			FY[j]=FY[j]-FK*(atomOfIPos.y-atomOfJPos.y)/r;
			FZ[i]=FZ[i]+FK*(atomOfIPos.z-atomOfJPos.z)/r;
			FZ[j]=FZ[j]-FK*(atomOfIPos.z-atomOfJPos.z)/r;
		}
	}

	for(i = 0;i < N; i++)
	{
		maxF = (FX[i]>maxF)?FX[i]:maxF;
		maxF = (FY[i]>maxF)?FY[i]:maxF;
		maxF = (FZ[i]>maxF)?FZ[i]:maxF;
	}

	free(PEN);
	return maxF;

}

PE_AtomParamter& FSPotentialEnergy::ReturnAtomParameter(ATOM_TYPE atom)
{
	return FS_AtomParamter(0,0,0,0,0,0,0);
}

PE_AtomParamter& FSPotentialEnergy::ReturnAtomParameter(ATOM_TYPE atom1,ATOM_TYPE atom2)
{
	return FS_AtomParamter(0,0,0,0,0,0,0);
}