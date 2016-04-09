#include "StdAfx.h"
#include "FileTool.h"
#include "InfoModel.h"

static string InfoFileName = "Info.txt";
static string ResultFileName = "Result.txt";
static string DiamondFileName = "Diamond.txt";
static string EnergyFileName = "Energy.txt";
static string AlgorithmFileName = "Algorithm.txt";

FileTool::FileTool(void)
{
}


FileTool::~FileTool(void)
{
}

void FileTool::CreatFolder( string &root )
{
	string subStr;
	int begin = 0;
	int end = 0;
	while ( end != root.length() )
	{
		end = root.find("//",begin);
		if ( end == string::npos)
			end = root.length();
		subStr = root.substr( 0, end);
		begin = end + 2;
		mkdir( subStr.c_str() );
	}
}

string FileTool::RootFolderName( string &baseRoot, PE_TYPE energy, Alloy &alloy, AlloyNum &num, int N )
{
	string fileName = baseRoot;
	RootAppend( fileName, PEnergyFolderName( energy ) );
	RootAppend( fileName, AlloyFolderName( alloy ) );
	RootAppend( fileName, SizeFolderName( N ) );
	RootAppend( fileName, AlloyNumFolderName( num ) );
	return fileName;
}

string FileTool::PEnergyFolderName( PE_TYPE type )
{
	string fileName = "";
	fileName = PE_EnergyTypeName( type );
	return fileName;
}

string FileTool::AlloyFolderName( Alloy &alloy )
{
	string fileName = "";
	for ( int i = 0; i < alloy.size(); i++ )
	{
		if ( i == 0 )
			fileName += AtomTypeName( alloy[0] );
		else
		{
			fileName += "-";
			fileName += AtomTypeName( alloy[i] );
		}
	}
	return fileName;
}

string FileTool::SizeFolderName( int N )
{
	char name[100];
	string fileName = "";

	sprintf( name, "%d", N );
	fileName += name;

	return fileName;
}

string FileTool::AlloyNumFolderName( AlloyNum &alloyNum )
{	
	char name[100];

	string fileName = "";
	for ( int i = 0; i < alloyNum.size(); i++ )
	{	
		sprintf( name, "%d", alloyNum[i] );
		if ( i != 0 )
			fileName += "-";
		fileName += name;
	}

	return fileName;
}

void FileTool::RootAppend( string &m, string n )
{
	m += "//";
	m += n;
}

string FileTool::FileNameWithType( string &root, FileType fileType )
{
	string fileName = root;

	switch( fileType )
	{
	case File_Info: { RootAppend( fileName, InfoFileName ); break; }
	case File_Result: { RootAppend( fileName, ResultFileName ); break; }
	case File_Diamond: { RootAppend( fileName, DiamondFileName ); break; }
	case File_Energy: { RootAppend( fileName, EnergyFileName ); break; }
	case File_Algorithm: { RootAppend( fileName, AlgorithmFileName ); break; }
	}

	return fileName;
}

void FileTool::SaveInfoFile( string &folderName, Clusters &cluster, PE_TYPE type )
{
	InfoModel model;
	InfoModel::ModelWithCluster( model, cluster );

	FileTool::CreatFolder( folderName );
	string infoFileName = FileTool::FileNameWithType( folderName, File_Info );

	model.SaveToFile( infoFileName );
}

void FileTool::SaveResultFile( string &folderName, Clusters &cluster )
{
	string resultFileName = FileTool::FileNameWithType( folderName, File_Result );

	FileTool::CreatFolder( folderName );

	ofstream file( resultFileName );
	int N = cluster.GetAtomsNumber();
	file<<N<<"\n";
	for ( int i = 0; i < N; i ++ )
	{
		Atom atom = cluster.GetAtomAtIndex( i );
		file<<setw(10)<<setiosflags(ios::left)<<atom.GetNote();
		file<<setw(20)<<setiosflags(ios::left)<<atom.GetPos().x;
		file<<setw(20)<<setiosflags(ios::left)<<atom.GetPos().y;
		file<<setw(20)<<setiosflags(ios::left)<<atom.GetPos().z<<endl;
	}

	file.close();
}

void FileTool::SaveDiamondFile( string &folderName, Clusters &cluster )
{
	string diamondFileName = FileTool::FileNameWithType( folderName, File_Diamond );

	FileTool::CreatFolder( folderName );

	ofstream file( diamondFileName );
	Alloy alloy = cluster.GetAlloy();
	int N = cluster.GetAtomsNumber();
	file<<N<<"\n";
	for ( int i = 0; i < N; i ++ )
	{
		Atom atom = cluster.GetAtomAtIndex( i );
		file<<setw(10)<<setiosflags(ios::left)<<AtomTypeName(alloy[atom.GetNote()]);
		file<<setw(20)<<setiosflags(ios::left)<<atom.GetPos().x;
		file<<setw(20)<<setiosflags(ios::left)<<atom.GetPos().y;
		file<<setw(20)<<setiosflags(ios::left)<<atom.GetPos().z<<endl;
	}

	file.close();
}