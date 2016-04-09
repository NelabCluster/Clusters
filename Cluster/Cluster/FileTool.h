#pragma once
#include "Common.h"

enum FileType
{
	File_Info,
	File_Result,
	File_Diamond,
	File_Energy,
	File_Algorithm
};

class FileTool
{
public:
	FileTool(void);
	~FileTool(void);

	static void CreatFolder( string &root );

	static string RootFolderName( string &baseRoot, PE_TYPE energyType, Alloy &alloy, AlloyNum &num, int N  );
	static string PEnergyFolderName( PE_TYPE type );
	static string AlloyFolderName( Alloy &alloy );
	static string SizeFolderName( int N );
	static string AlloyNumFolderName( AlloyNum &alloyNum );

	static string FileNameWithType( string &root, FileType fileType );
	
	static void SaveInfoFile( string &folderName, Clusters &cluster, PE_TYPE type );
	static void SaveResultFile( string &folderName, Clusters &cluster );
	static void SaveDiamondFile( string &folderName, Clusters &cluster );

	static void RootAppend( string &m, string n );
};



