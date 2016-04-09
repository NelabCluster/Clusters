#include "StdAfx.h"
#include "InitTool.h"
#include "BaseTool.h"
#include "LocalTool.h"

void InitTool::RandInCubic(Clusters& cluster, double length)
{
	for ( int i = 0; i < cluster.GetAtomsNumber(); i++ )
	{	
		double x = (RANDI-0.5) * length;
		double y = (RANDI-0.5) * length;
		double z = (RANDI-0.5) * length;
		cluster.SetPosOfAtomIndex(i,AtomPos(x,y,z));
	}
}

void InitTool::RandInSphere(Clusters& cluster, double radius)
{
	for ( int i = 0; i < cluster.GetAtomsNumber(); i++)
	{
		double r = RANDI * radius;
		double theta = RANDI * M_PI;
		double phi = RANDI * 2 * M_PI;
		cluster.SetPosOfAtomIndex(i,AtomPos(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)));
	}
}

void InitTool::GrowFromSmallerCluster(Clusters& cluster, Clusters& smallClusters)
{
	assert(cluster.GetAtomsNumber() > smallClusters.GetAtomsNumber());

	int diff = cluster.GetAtomsNumber() - smallClusters.GetAtomsNumber();

	smallClusters.MoveCenterToOrigin();

	for (int i = 0; i < smallClusters.GetAtomsNumber(); i++)
	{
		Atom atom = smallClusters.GetAtomAtIndex(i);
		cluster.SetAtomOfAtomIndex(i,atom);
	}

	double maxDisToOrigin = DBL_MIN;
	for (int i = 0; i < smallClusters.GetAtomsNumber(); i++)
	{
		Atom atom = smallClusters.GetAtomAtIndex(i);
		AtomPos pos = atom.GetPos();
		double disToOrigin = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
		
		if (disToOrigin > maxDisToOrigin) maxDisToOrigin = disToOrigin;
	}
	maxDisToOrigin = sqrt(maxDisToOrigin);

	for (int i = 0; i < diff; i++)
	{
		double r = RANDI * maxDisToOrigin;
		double theta = RANDI * M_PI;
		double phi = RANDI * M_PI * 2;
		AtomPos pos(XFROMPOLE(r,theta,phi),YFROMPOLE(r,theta,phi),ZFROMPOLE(r,theta,phi));
		cluster.SetPosOfAtomIndex(i,pos);
	}
}

void InitTool::RandType(Clusters& cluster)
{
	int N = cluster.GetAtomsNumber();
	AlloyNum alloyNum = cluster.GetAlloyNum();
	int atomTypeCount = alloyNum.size();

	vector<int> sumNumberBeforeAtom(atomTypeCount,0);
	sumNumberBeforeAtom[0] = alloyNum[0];
	for(int i=1;i<atomTypeCount;i++)
		sumNumberBeforeAtom[i] = sumNumberBeforeAtom[i-1] + alloyNum[i];

	vector<int> note = BaseTool::RandPerm(N,N);

	for(int i=0;i<N;i++){
		int index = 0;
		while( note[i] >= sumNumberBeforeAtom[index] )
			index ++;
		note[i] = index;
	}

	for (int i=0; i < N; i++)
		cluster.SetNoteOfAtomIndex(i,note[i]);
}

void InitTool::AtomTypeSwap( Clusters &cluster, Clusters &neighCluster, PE_TYPE type )
{
	if ( cluster.GetAtomsNumber() != neighCluster.GetAtomsNumber() )
		return;

	if ( cluster.GetAlloy().size() != 2 || neighCluster.GetAlloy().size() != 2 )
		return;

	int N = cluster.GetAtomsNumber();

	int diff = cluster.GetAlloyNum()[0] - neighCluster.GetAlloyNum()[0];

	if ( abs(diff) != 1 )
		return;

	int note = (diff == 1)?0:1;

	Clusters temp = cluster;
	cluster.SetEnergy( 0 );
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++)
			temp.SetAtomOfAtomIndex( j, neighCluster.GetAtomAtIndex(j) );
		Atom atom = temp.GetAtomAtIndex( i );
		if ( atom.GetNote() == (1-note) )
		{
			temp.SetNoteOfAtomIndex( i, note );
			LocalTool::LocalMinimizeClusters( temp, type );
			if ( temp.GetEnergy() < cluster.GetEnergy() )
				cluster = temp;
		}
	}
}
