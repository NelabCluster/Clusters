#include "StdAfx.h"
#include "DEBaseClusters.h"


DEBaseClusters::DEBaseClusters(int N,vector<ATOM_TYPE> types,vector<int> numbers)
	:Clusters(N,types,numbers)
{
}

DEBaseClusters::~DEBaseClusters(void)
{
}

DEBaseClusters operator+(DEBaseClusters& left,DEBaseClusters& right)
{
	DEBaseClusters newCluster(left);
	for (int i = 0; i < newCluster.GetAtomsNumber(); i ++)
	{
		AtomPos atomRightPos = right.GetAtomAtIndex(i).GetPos();
		AtomPos atomLeftPos = left.GetAtomAtIndex(i).GetPos();
		double x = atomLeftPos.x + atomRightPos.x;
		double y = atomLeftPos.y + atomRightPos.y;
		double z = atomLeftPos.z + atomRightPos.z;
		newCluster.SetPosOfAtomIndex(i,AtomPos(x,y,z));
	}
	return newCluster;
}

DEBaseClusters operator-(DEBaseClusters& left,DEBaseClusters& right)
{
	DEBaseClusters newCluster(left);
	for (int i = 0; i < newCluster.GetAtomsNumber(); i ++)
	{
		AtomPos atomRightPos = right.GetAtomAtIndex(i).GetPos();
		AtomPos atomLeftPos = left.GetAtomAtIndex(i).GetPos();
		double x = atomLeftPos.x - atomRightPos.x;
		double y = atomLeftPos.y - atomRightPos.y;
		double z = atomLeftPos.z - atomRightPos.z;
		newCluster.SetPosOfAtomIndex(i,AtomPos(x,y,z));
	}
	return newCluster;
}

DEBaseClusters operator*(double value,DEBaseClusters& right)
{
	DEBaseClusters newCluster(right);
	for (int i = 0; i < newCluster.GetAtomsNumber(); i ++)
	{
		AtomPos atomRightPos = right.GetAtomAtIndex(i).GetPos();
		double x = value * atomRightPos.x;
		double y = value * atomRightPos.y;
		double z = value * atomRightPos.z;
		newCluster.SetPosOfAtomIndex(i,AtomPos(x,y,z));
	}
	return newCluster;
}

DEBaseClusters operator*(DEBaseClusters& left,double value)
{
	return value * left;
}