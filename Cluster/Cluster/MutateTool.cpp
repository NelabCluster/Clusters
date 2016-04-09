#include "StdAfx.h"
#include "MutateTool.h"

void MutateTool::AtomDisplacement(Clusters& cluster,int number)
{
	int N = cluster.GetAtomsNumber();
	vector<bool> flag;
	flag.resize(N,false);

	cluster.MoveCenterToOrigin();

	while( number != 0 )
	{
		int index = RANDIUINT(N);
		if ( flag[index] == false )
		{
			flag[index] = true;
			number --;
		}
	}

	double maxR = 0;
	for ( int i = 0; i < N; i++ )
	{
		AtomPos pos = cluster.GetAtomAtIndex( i ).GetPos();
		double r = sqrt( pos.x * pos.x + pos.y * pos.y + pos.z * pos.z );
		if ( r > maxR)
			maxR = r;
	}

	//double r0 = AtomPara::NearestNeighborSeparation(cluster.GetAlloy());
	double length = maxR*2;//r0 * pow( N * 1.0, 1.0/3 );
	for (int i = 0; i < N; i++ )
	{
		if (flag[i])
		{
			double x = (RANDI-0.5) * length;
			double y = (RANDI-0.5) * length;
			double z = (RANDI-0.5) * length;
			cluster.SetPosOfAtomIndex(i,AtomPos(x,y,z));
		}
	}
}

void MutateTool::Twisting(Clusters& cluster)
{
	int N = cluster.GetAtomsNumber();
	vector<bool> flag;
	flag.resize(N,false);

	cluster.MoveCenterToOrigin();

	for (int i = 0; i < N; i++ )
	{
		AtomPos pos = cluster.GetAtomAtIndex( i ).GetPos();
		if ( pos.z > 0 )
			flag[i] = true;
	}
	cluster.RotateWithVectorAndAngle(flag,0,0,1,RANDI*M_PI*2);
}

void MutateTool::Replacement(Clusters& cluster)
{
	int N = cluster.GetAtomsNumber();
	cluster.MoveCenterToOrigin();



	double r0 = AtomPara::NearestNeighborSeparation(cluster.GetAlloy());

	double length = r0 * pow( N * 1.0, 1.0/3 );
	for ( int i = 0; i < N; i++ )
	{	
		double x = (RANDI-0.5) * length;
		double y = (RANDI-0.5) * length;
		double z = (RANDI-0.5) * length;
		cluster.SetPosOfAtomIndex(i,AtomPos(x,y,z));
	}
}

void MutateTool::AtomRandSwap(Clusters& cluster, int times)
{
	int N = cluster.GetAtomsNumber();

	while(times!=0)
	{
		int i = RANDIUINT(N);
		int j = RANDIUINT(N);
		
		while(cluster.GetAtomAtIndex(i).GetNote() == cluster.GetAtomAtIndex(j).GetNote())
		{
			i = RANDIUINT(N);
			j = RANDIUINT(N);
		}

		int t1 = cluster.GetAtomAtIndex(i).GetNote();
		int t2 = cluster.GetAtomAtIndex(j).GetNote();
		cluster.SetNoteOfAtomIndex(j,t1);
		cluster.SetNoteOfAtomIndex(i,t2);

		times --;
	}
}