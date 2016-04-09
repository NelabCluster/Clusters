// Cluster.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Common.h"
#include "Tool.h"
#include "BCGA.h"
#include "AlgorithmManage.h"


//-3.7513
void TestEnergy()
{
	Clusters cluster(3,ATOM_Zn);
	cluster.SetAtomOfAtomIndex(0,Atom(0,AtomPos(-0.0025257515,-1.1132616080,-0.1490790579)));
	cluster.SetAtomOfAtomIndex(1,Atom(0,AtomPos(0.6823000831,-0.6720141616,0.1308423179)));
	cluster.SetAtomOfAtomIndex(2,Atom(0,AtomPos(0.2779619653,-1.2199766042,0.6583738879)));
	PotentialEnergy * pe = PotentialEnergy::initWithType(PE_Gupta);
	double E = pe->EnergyValue(cluster);
	cout<<E<<endl; 
	delete pe;
}

void TestLocal()
{
	Clusters cluster(3,ATOM_Zn);
	cluster.SetAtomOfAtomIndex(0,Atom(0,AtomPos(1,0,0)));
	cluster.SetAtomOfAtomIndex(1,Atom(0,AtomPos(-1,0,0)));
	cluster.SetAtomOfAtomIndex(2,Atom(0,AtomPos(0,1,0)));
	//InitTool::RandInCubic(cluster,2);
	PotentialEnergy * pe = PotentialEnergy::initWithType(PE_Gupta);
	double E = pe->EnergyValue(cluster);
	double F = pe->ForceValue(cluster);
	cout<<E<<endl;
	cout<<F<<endl;
	E = LocalTool::LocalMinimizeClusters(cluster,PE_Gupta);
	cout<<E<<endl; 
	delete pe;
}

void TestFileName()
{
	Alloy alloy;
	alloy.push_back( ATOM_Cu );
	alloy.push_back( ATOM_Au );
	string alloyFileName = FileTool::AlloyFolderName( alloy );
	cout<<alloyFileName<<endl;

	AlloyNum alloyNum;
	alloyNum.push_back(21);
	alloyNum.push_back(7);
	string alloyNumFileName = FileTool::AlloyNumFolderName( alloyNum );
	cout<<alloyNumFileName<<endl;

	PE_TYPE energy = PE_Gupta;

	//FileTool::CreatRootDocument( energy, alloy, alloyNum, 28 );
	(ClusterPool::Share())->GetEnergy( energy, alloy, alloyNum, 28 );

}

void TestGetCluster()
{
	Clusters one( 28, ATOM_Cu, ATOM_Au, 21 );
	PE_TYPE type = PE_Gupta;

	ClusterPool::Share()->SetCluster( one, type );

	Clusters two = one;
	ClusterPool::Share()->GetCluster( two, type );

	ClusterPool::Share()->SetDiamomd( one, type );
	
	cout<<"Test end"<<endl;
}

void TestRandINT()
{
	int all = 0, m = 0;

	for ( int i = 0; i < 1000; i++ )
	{
		int a = RANDIUINT(3);
		if ( a == 0 )
			m ++;
		all ++;
	}
	cout<<m<<"\t"<<all<<endl;
}

void TestFingerprint()
{
	Clusters cluster(10,ATOM_Cu);

	int N = cluster.GetAtomsNumber();
	double r0 = AtomPara::NearestNeighborSeparation( cluster.GetAlloy() );
	double radius = r0 * pow( N * 1.0, 1.0/3 );
	InitTool::RandInSphere( cluster, radius );

	LocalTool::LocalMinimizeClusters( cluster, PE_Gupta );
	PropertyTool::FingerprintFunction( cluster );
}

void PreExperiment()
{
	int N = 28,m = 0;

	for ( int m = 0; m <= N; m++)
	{
		string root = "Pre";
		BCGAParameters para(N,PE_Gupta,ATOM_Cu,ATOM_Au,m,root);
		para.gen = (N<20)?100:150;
		para.popSize = (N<20)?20:30;
		if ( m == 0 || m == N )
			para.mutate = BCGA_Mutate_Displace;
		BCGAAlgorithm alg(para);
		AlgorithmManage manage( &alg );
		manage.start(10);
	}

	for ( int m = N; m >=0; m--)
	{
		string root = "Pre";
		BCGAParameters para(N,PE_Gupta,ATOM_Cu,ATOM_Au,m,root);
		para.gen = (N<20)?100:150;
		para.popSize = (N<20)?20:30;
		if ( m == 0 || m == N )
			para.mutate = BCGA_Mutate_Displace;
		BCGAAlgorithm alg(para);
		AlgorithmManage manage( &alg );
		manage.start(10);
	}
}

void PopsizeExperiment()
{
	int N = 16,m = 4,popSize=30;
	cin>>N>>popSize;
	int o = N / 4;
	//cin>>N>>m>>popSize;
	//if ( N == 15 ) m = 1;
	//if ( N == 16 ) m = 4;
	//if ( N == 28 ) m = 9;
	
	for ( int i = 1; i <= 3; i++)
	{
		m = o * i;
		string root = "POPSIZE";
		FileTool::RootAppend(root, BaseTool::IntToString(N) );
		FileTool::RootAppend(root, BaseTool::IntToString(popSize));
		FileTool::RootAppend(root,BaseTool::IntToString( m ) );
		for ( int n = 0; n < 10; n++ )
		{
			string tempRoot = root;
			FileTool::RootAppend( tempRoot, BaseTool::IntToString( n ));
			BCGAParameters para(N,PE_Gupta,ATOM_Cu,ATOM_Au,m,tempRoot);
			para.init = BCGA_Init_Rand;

			para.gen = (N<20)?100:150;

			para.popSize = popSize;
		//	para.crossoverRate = 0;
//			para.muteRate = 1.0;
			if ( m == 0 || m == N )
				para.mutate = BCGA_Mutate_Displace;

			BCGAAlgorithm alg(para);
			AlgorithmManage manage( &alg );
			manage.start();
		}
	
	}
}

void InitTypeExperiment()
{
	
}

void MatRateExperiment()
{
	double crossoverRates[] = {0.0,0.2,0.4,0.6,0.8,1.0};

	int N = 28,m = 4;
	int o = N / 4;

	for ( int i = 1; i <= 3; i++)
	{
		m = o * i;
		for ( int j = 0; j < 6; j ++)
		{
			double crossoverRate = crossoverRates[j];
			string root = "MateRate";
			FileTool::RootAppend(root, BaseTool::IntToString(N) );
			FileTool::RootAppend(root, BaseTool::DoubleToString(crossoverRate,1));
			FileTool::RootAppend(root,BaseTool::IntToString( m ) );
			for ( int n = 0; n < 10; n++ )
			{
				string tempRoot = root;
				FileTool::RootAppend( tempRoot, BaseTool::IntToString( n ));
				BCGAParameters para(N,PE_Gupta,ATOM_Cu,ATOM_Au,m,tempRoot);
				para.gen = (N<20)?100:150;
				para.popSize = 20;
				para.init = BCGA_Init_Rand;
				para.crossoverRate = crossoverRate;
				if ( m == 0 || m == N )
					para.mutate = BCGA_Mutate_Displace;
				BCGAAlgorithm alg(para);
				AlgorithmManage manage( &alg );
				manage.start();
			}
		}
	}
}

void MateTypeExperiment()
{
	int N = 28,m = 4;
	int o = N / 4;

	for ( int i = 1; i <= 3; i++)
	{
		m = o * i;
		for ( int j = 0; j < 2; j ++)
		{
			string mateTypeName = (j==0)?"Value":"CutAndSplice"; 
			BCGA_Crossover mateType = (j==0)?BCGA_Crossover_Value:BCGA_Crossover_PlaneCutSplice;

			string root = "MateType";
			FileTool::RootAppend(root, BaseTool::IntToString(N) );
			FileTool::RootAppend(root, mateTypeName);
			FileTool::RootAppend(root,BaseTool::IntToString( m ) );
			for ( int n = 0; n < 10; n++ )
			{
				string tempRoot = root;
				FileTool::RootAppend( tempRoot, BaseTool::IntToString( n ));
				BCGAParameters para(N,PE_Gupta,ATOM_Cu,ATOM_Au,m,tempRoot);
				para.gen = (N<20)?100:150;
				para.init = BCGA_Init_Rand;
				para.crossover = mateType;
				if ( m == 0 || m == N )
					para.mutate = BCGA_Mutate_Displace;
				BCGAAlgorithm alg(para);
				AlgorithmManage manage( &alg );
				manage.start();
			}
		}
	}
}

void MutateTypeExperiment()
{
	int N = 28,m = 4;
	int o = N / 4;

	for ( int i = 1; i <= 3; i++)
	{
		m = o * i;
		for ( int j = 0; j < 6; j ++)
		{
			string name; 
			int type;

			switch (j)
			{
			case 0: name="displace"; type = BCGA_Mutate_Displace; break;
			case 1: name="twitst"; type = BCGA_Mutate_Twitst; break;
			case 2: name="replace"; type = BCGA_Mutate_Replace; break;
			case 3: name="swap"; type = BCGA_Mutate_Swap; break;
			case 4: name="displaceswap"; type = BCGA_Mutate_Displace | BCGA_Mutate_Swap; break;
			case 5: name="all"; type = BCGA_Mutate_Displace | BCGA_Mutate_Twitst | BCGA_Mutate_Replace | BCGA_Mutate_Swap; break;
			}

			string root = "MutateType";
			FileTool::RootAppend(root, BaseTool::IntToString(N) );
			FileTool::RootAppend(root, name);
			FileTool::RootAppend(root,BaseTool::IntToString( m ) );
			for ( int n = 0; n < 10; n++ )
			{
				string tempRoot = root;
				FileTool::RootAppend( tempRoot, BaseTool::IntToString( n ));
				BCGAParameters para(N,PE_Gupta,ATOM_Cu,ATOM_Au,m,tempRoot);
				para.gen = (N<20)?100:150;
				para.init = BCGA_Init_Rand;
				para.mutate = type;
				BCGAAlgorithm alg(para);
				AlgorithmManage manage( &alg );
				manage.start();
			}
		}
	}
}

void MutateRateExperiment()
{
	double mutateRates[] = {0.1,0.3,0.5,0.7,0.9};

	int N = 28,m = 4;
	int o = N / 4;

	for ( int i = 1; i <= 3; i++)
	{
		m = o * i;
		for ( int j = 0; j < 5; j ++)
		{
			double mutateRate = mutateRates[j];
			string root = "MutateRate";
			FileTool::RootAppend(root, BaseTool::IntToString(N) );
			FileTool::RootAppend(root, BaseTool::DoubleToString(mutateRate,1));
			FileTool::RootAppend(root,BaseTool::IntToString( m ) );
			for ( int n = 0; n < 10; n++ )
			{
				string tempRoot = root;
				FileTool::RootAppend( tempRoot, BaseTool::IntToString( n ));
				BCGAParameters para(N,PE_Gupta,ATOM_Cu,ATOM_Au,m,tempRoot);
				para.gen = (N<20)?100:150;
				para.init = BCGA_Init_Rand;
				para.muteRate = mutateRate;
				if ( m == 0 || m == N )
					para.mutate = BCGA_Mutate_Displace;
				BCGAAlgorithm alg(para);
				AlgorithmManage manage( &alg );
				manage.start();
			}
		}
	}
}

void MutateNumberExperiment()
{
	int N = 28,m = 4;
	int o = N / 4;

	for ( int i = 1; i <= 3; i++)
	{
		m = o * i;
		for ( int j = 0; j < 2; j ++)
		{
			string name; 
			int type;

			switch (j)
			{
			case 0: name="displace"; type = BCGA_Mutate_Displace; break;
			case 1: name="swap"; type = BCGA_Mutate_Swap; break;
			}

			for ( int num = 1; num <=4; num ++ )
			{
				string root = "MutateNumber";
				FileTool::RootAppend(root, BaseTool::IntToString(N) );
				FileTool::RootAppend(root, name);
				FileTool::RootAppend(root,BaseTool::IntToString( m ) );
				FileTool::RootAppend(root, BaseTool::IntToString( num ));
				for ( int n = 0; n < 10; n++ )
				{
					string tempRoot = root;
					FileTool::RootAppend( tempRoot, BaseTool::IntToString( n ));
					BCGAParameters para(N,PE_Gupta,ATOM_Cu,ATOM_Au,m,tempRoot);
					para.gen = (N<20)?100:150;
					para.init = BCGA_Init_Rand;
					para.mutate = type;
					para.mutateNum = num;
						BCGAAlgorithm alg(para);
					AlgorithmManage manage( &alg );
					manage.start();
				}
			}
		}
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	srand((unsigned)time(NULL));

	#pragma region 测试函数
		//TestEnergy();
		//TestLocal();
		//TestFileName();
		//TestGetCluster();
		//TestRandINT();
		//TestFingerprint();
	#pragma endregion 测试函数

	srand((unsigned)time(NULL));

	PreExperiment();

	//PopsizeExperiment();
	//InitTypeExperiment();

	//MatRateExperiment();
	//MateTypeExperiment();

	//MutateTypeExperiment();
	//MutateRateExperiment();
	//MutateNumberExperiment();

	system("pause");
	return 0;
}
