#include "StdAfx.h"
#include "AlgorithmManage.h"
#include "FileTool.h"
#include "BaseTool.h"

AlgorithmManage::AlgorithmManage(void) { }

AlgorithmManage::~AlgorithmManage(void) { }

AlgorithmManage::AlgorithmManage( IAlgorithm* a ):_algorithm(a) { }

void AlgorithmManage::start()
{	
	this->_Start( );
}

void AlgorithmManage::start(int repeat)
{	
	int success = 0;
	int times = repeat;

	while(repeat!=0)
	{
		//string folderName = _root;
		//string secondFolderName = BaseTool::IntToString( times-repeat );
		//FileTool::RootAppend( folderName, secondFolderName );

		this->_Start( );
		if ( fabs(_algorithm->GetBestEnergy() - (-126.592)) < 0.001 )
		{
			success ++;
		}
		repeat -- ;
		cout<<success<<"/"<<times-repeat<<endl;
	}

}

void AlgorithmManage::_Start( )
{
	srand((unsigned)time(NULL));
	_algorithm->Initialization();

	while (!_algorithm->EndCondition())
	{
		_algorithm->Process();
		_algorithm->PrintProgress();
		_algorithm->SaveProgress( );
	}
	_algorithm->Ending();
}

