#include "StdAfx.h"
#include "ClusterPool.h"
#include "InfoModel.h"

static string PoolFileName = "Clusters";

void ClusterPool::SaveCluster( Clusters &cluster, PE_TYPE energyType )
{
	SetInfo( cluster, energyType );
	SetCluster( cluster, energyType );
	SetDiamomd( cluster, energyType );
}


double ClusterPool::GetEnergy( PE_TYPE energyType, Alloy &alloy, AlloyNum &num, int N )
{
	string folderName = FileTool::RootFolderName( PoolFileName, energyType, alloy, num, N );
	string infoFileName = FileTool::FileNameWithType( folderName, File_Info );

	InfoModel model;
	
	if( !InfoModel::ModelFromFile( model, infoFileName ) )
		return 0;
	
	double energy = model.energy;

	return energy;
}

void ClusterPool::SetInfo( Clusters &cluster, PE_TYPE energyType )
{
	string folderName = FileTool::RootFolderName( PoolFileName, energyType, cluster.GetAlloy(), cluster.GetAlloyNum(), cluster.GetAtomsNumber() );
	FileTool::SaveInfoFile( folderName, cluster, energyType );
}

bool ClusterPool::GetCluster( Clusters &cluster, PE_TYPE energyType )
{
	string folderName = FileTool::RootFolderName( PoolFileName, energyType, cluster.GetAlloy(), cluster.GetAlloyNum(), cluster.GetAtomsNumber() );
	string resultFileName = FileTool::FileNameWithType( folderName, File_Result );

	ifstream file( resultFileName );

	if ( !file.is_open() )
		return false;

	int N;
	file>>N;

	for ( int i = 0; i < N; i++ )
	{
		int note;
		double x,y,z;
		file>>note>>x>>y>>z;
		Atom atom( note, AtomPos(x,y,z) );
		cluster.SetAtomOfAtomIndex( i, atom );
	}

	file.close();
}

void ClusterPool::SetCluster( Clusters &cluster,  PE_TYPE energyType  )
{
	string folderName = FileTool::RootFolderName( PoolFileName, energyType, cluster.GetAlloy(), cluster.GetAlloyNum(), cluster.GetAtomsNumber() );
	FileTool::SaveResultFile( folderName, cluster );
}

void ClusterPool::SetDiamomd( Clusters &cluster,  PE_TYPE energyType  )
{
	string folderName = FileTool::RootFolderName( PoolFileName, energyType, cluster.GetAlloy(), cluster.GetAlloyNum(), cluster.GetAtomsNumber() );
	FileTool::SaveDiamondFile( folderName, cluster );
}
