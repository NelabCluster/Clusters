#include "StdAfx.h"
#include "InfoModel.h"
#include "FileTool.h"

InfoModel::InfoModel(void)
{
}


InfoModel::~InfoModel(void)
{
}

void InfoModel::ModelWithCluster( InfoModel &model, Clusters &cluster )
{
	model.energy = cluster.GetEnergy();
	model.N = cluster.GetAtomsNumber();
	model.alloy = cluster.GetAlloy();
	model.alloyNum = cluster.GetAlloyNum();
}

bool InfoModel::ModelFromFile( InfoModel& model, string &fileName )
{
	ifstream file( fileName );

	if ( !file.is_open())
		return false;

	file>>model.energy;

	file.close();

	return true;
}

void InfoModel::SaveToFile( string &fileName )
{
	ofstream file( fileName );

	file<<energy<<endl;

	file<<N<<endl;

	file<<FileTool::AlloyFolderName( alloy )<<endl;

	file<<FileTool::AlloyNumFolderName( alloyNum )<<endl;

	file.close();
}
