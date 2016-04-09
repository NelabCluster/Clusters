#include "StdAfx.h"
#include "CrossOverTool.h"
#include "PotentialEnergy.h"

void CrossOverTool::PlaneCutSplice(Clusters& mother, Clusters& father)
{
	Clusters& offSpring1Clusters = mother;
	Clusters& offSpring2Clusters = father;	

	offSpring1Clusters.RotateOriginVector(RANDI * M_PI, RANDI * 2 * M_PI, RANDI * 2 * M_PI);
	offSpring2Clusters.RotateOriginVector(RANDI * M_PI, RANDI * 2 * M_PI, RANDI * 2 * M_PI);

	offSpring1Clusters.AtomsOrderByZ();
	offSpring2Clusters.AtomsOrderByZ();

	int randNumber = RANDIUINT(offSpring1Clusters.GetAtomsNumber()-1);
	//cout<<randNumber<<endl;

	offSpring1Clusters.MoveCenterToOrigin();
	offSpring2Clusters.MoveCenterToOrigin();

	//cout<<"........"<<endl;
	//for (int i =0; i < offSpring1Clusters.GetAtomsNumber(); i++)
	//{
	//	cout<<offSpring1Clusters.GetAtomAtIndex(i).GetNote()<<"\t"<<offSpring1Clusters.GetAtomAtIndex(i).GetPos().z;
	//	cout<<"\t";
	//	cout<<offSpring2Clusters.GetAtomAtIndex(i).GetNote()<<"\t"<<offSpring2Clusters.GetAtomAtIndex(i).GetPos().z<<endl;
	//}

	double dis1 = (offSpring1Clusters.GetAtomAtIndex(randNumber).GetPos().z + offSpring1Clusters.GetAtomAtIndex(randNumber+1).GetPos().z) / 2;
	offSpring1Clusters.MoveCenter(0,0,-dis1);
	double dis2 = (offSpring2Clusters.GetAtomAtIndex(randNumber).GetPos().z + offSpring2Clusters.GetAtomAtIndex(randNumber+1).GetPos().z) / 2;
	offSpring2Clusters.MoveCenter(0,0,-dis2);

	//cout<<"........"<<endl;
	//for (int i =0; i < offSpring1Clusters.GetAtomsNumber(); i++)
	//{
	//	cout<<offSpring1Clusters.GetAtomAtIndex(i).GetNote()<<"\t"<<offSpring1Clusters.GetAtomAtIndex(i).GetPos().z;
	//	cout<<"\t";
	//	cout<<offSpring2Clusters.GetAtomAtIndex(i).GetNote()<<"\t"<<offSpring2Clusters.GetAtomAtIndex(i).GetPos().z<<endl;
	//}

	Clusters temp0(offSpring1Clusters);
	PlaneCutSpliceOperate(offSpring1Clusters,offSpring2Clusters,randNumber);
	PlaneCutSpliceOperate(offSpring2Clusters,temp0,randNumber);

	//cout<<"........"<<endl;
	//for (int i =0; i < offSpring1Clusters.GetAtomsNumber(); i++)
	//{
	//	cout<<offSpring1Clusters.GetAtomAtIndex(i).GetNote()<<"\t"<<offSpring1Clusters.GetAtomAtIndex(i).GetPos().z;
	//	cout<<"\t";
	//	cout<<offSpring2Clusters.GetAtomAtIndex(i).GetNote()<<"\t"<<offSpring2Clusters.GetAtomAtIndex(i).GetPos().z<<endl;
	//}
}

void CrossOverTool::PlaneCutSpliceOperate(Clusters& top,Clusters& bottom,int cutPoint)
{
	AlloyNum alloyNum = top.GetAlloyNum();
	int typeNumber = alloyNum.size();

	for (int i = 0; i <= cutPoint; i++)
	{
		Atom atom = top.GetAtomAtIndex(i);
		alloyNum[atom.GetNote()] --;
	}

	int index0 = cutPoint+1;
	for (int i = 0; i < typeNumber; i++)
	{
		int index = bottom.GetAtomsNumber() - 1;
		while ( alloyNum[i] != 0 )
		{
			Atom atom = bottom.GetAtomAtIndex(index);
			if ( atom.GetNote() == i ) 
			{
				top.SetAtomOfAtomIndex(index0,atom);
				alloyNum[i] --;
				index0 ++;
			}
			index --;
		}
	}
}