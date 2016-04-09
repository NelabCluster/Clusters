#include "StdAfx.h"
#include "BCGAAlgorithm.h"
#include "Common.h"
#include "Tool.h"


#pragma region 构造函数

BCGAAlgorithm::BCGAAlgorithm(BCGAParameters& para)
	:_para(para),_iteration(0)
{
	_pe = PotentialEnergy::initWithType(_para.energyType);
	for (int i = 0; i < _para.popSize; i++)
	{
		BCGAClusters cluster(_para.N,_para.atomTypes,_para.atomNumber);
		_pop.push_back(cluster);
	}
}

BCGAAlgorithm::~BCGAAlgorithm(void)
{
	delete _pe;
}

#pragma endregion 构造函数

#pragma region IAlgorithm 协议

void BCGAAlgorithm::Initialization()
{
	_iteration = 0;
	_convergenceCount = 0;
	_lastPopBestEnergy = 0;
	_beginTime = clock();

	switch ( _para.init )
	{
	case BCGA_Init_Rand:
		this->RandInitialization();
		break;
	case BCGA_Init_Grow:
		this->GrowInitialization();
		break;
	case BCGA_Init_Swap:
		this->SwapInitialization();
		break;
	default:
		this->RandInitialization();
	}

	sort(_pop.begin(),_pop.end(),SortClusterByEnergy);

	this->SaveProgress();
}

void BCGAAlgorithm::Process()
{
	_iteration ++;

	//this->PrintfEnergy(_pop);
	vector<BCGAClusters> selected = this->Selection();
	//	cout<<"select"<<endl; this->PrintfEnergy(selected);

	vector<BCGAClusters> mat = this->CrossOver(selected);
	//vector<BCGAClusters> mat;
	//cout<<"mat"<<endl; this->PrintfEnergy(mat);

	vector<BCGAClusters> mut = this->Mutation(_pop);
	//cout<<"mut"<<endl; this->PrintfEnergy(mut);

	this->Choose1(mat,mut);
	//	this->PrintfEnergy(_pop);

}

bool BCGAAlgorithm::EndCondition()
{
	//double bestEnergy = (_pop.end()-1)->GetEnergy();
	//double worstEnergy = _pop.begin()->GetEnergy();

	//if ( bestEnergy == worstEnergy )
	//	_convergenceCount ++;
	//else
	//	_convergenceCount = 0;

	//return (_convergenceCount>30);

	return !(_iteration < _para.gen);
}

double BCGAAlgorithm::GetBestEnergy()
{
	return _pop[_pop.size()-1].GetEnergy();
}

void BCGAAlgorithm::PrintProgress()
{
	double aveEnergy = 0;
	for (int i = 0; i < _pop.size(); i++ )
		aveEnergy += _pop[i].GetEnergy();
	aveEnergy /= _pop.size();

	cout<<_iteration;
	cout<<_pop[_pop.size()-1].GetEnergy() ;
	cout<<_pop[0].GetEnergy() ;
	cout<<aveEnergy<<endl;
}

void BCGAAlgorithm::SaveProgress()
{
	this->SaveBestCluster();

	Clusters &best = _pop[_pop.size()-1];
	if ( best.GetEnergy() < _lastPopBestEnergy )
	{
		_lastPopBestEnergy = best.GetEnergy();
		FileTool::SaveInfoFile( _para.rootFolderName, best, _para.energyType );
		FileTool::SaveResultFile( _para.rootFolderName, best );
		FileTool::SaveDiamondFile( _para.rootFolderName, best );
	}
	this->SaveEnergyFile();
}

void BCGAAlgorithm::Ending()
{
	_endTime = clock();
	clock_t diff = _endTime - _beginTime;

	string folderName = _para.rootFolderName;
	FileTool::CreatFolder( folderName );
	string algorithmFileName = FileTool::FileNameWithType( folderName, File_Algorithm );

	ofstream file( algorithmFileName, ios::app );
	file<<diff<<endl;
	file.close();

}

#pragma endregion IAlgorithm 协议

#pragma region 不同的初始化

void BCGAAlgorithm::RandInitialization()
{
	double r0 = AtomPara::NearestNeighborSeparation(_para.atomTypes);
	double length = r0 * pow(_para.N * 1.0, 1.0/3);
	for (int i = 0; i < _para.popSize; i++)
	{
		BCGAClusters& cluster = _pop[i];
		InitTool::RandInCubic(cluster,length);
		InitTool::RandType(cluster);
		this->Evaluate(cluster);
	}
}

void BCGAAlgorithm::GrowInitialization()
{
	int m = 1;
	AlloyNum smallAlloyNum = _para.atomNumber;
	smallAlloyNum[0] --;
	Clusters smallCluster( _para.N - 1, _para.atomTypes, smallAlloyNum );
	if  ( ClusterPool::Share()->GetCluster( smallCluster, _para.energyType ) )
	{
		InitTool::GrowFromSmallerCluster( _pop[0], smallCluster );
		this->Evaluate( _pop[0] );
	}
	else
		m = 0;

	double r0 = AtomPara::NearestNeighborSeparation(_para.atomTypes);
	double length = r0 * pow(_para.N * 1.0, 1.0/3);
	for (int i = m; i < _para.popSize; i++)
	{
		BCGAClusters& cluster = _pop[i];
		InitTool::RandInCubic(cluster,length);
		InitTool::RandType(cluster);
		this->Evaluate(cluster);
	}
}

void BCGAAlgorithm::SwapInitialization()
{
	int m = 0;

	{
		AlloyNum neighAlloyNum = _para.atomNumber;
		neighAlloyNum[0]--;
		neighAlloyNum[1]++;
		Clusters neighCluster( _para.N, _para.atomTypes, neighAlloyNum );
		if ( ClusterPool::Share()->GetCluster( neighCluster, _para.energyType ) )
		{
			Clusters &one = _pop[m];
			InitTool::AtomTypeSwap( one, neighCluster, _para.energyType );
			m ++;
		} 
	}

	{
		AlloyNum neighAlloyNum = _para.atomNumber;
		neighAlloyNum[0]++;
		neighAlloyNum[1]--;
		Clusters neighCluster( _para.N, _para.atomTypes, neighAlloyNum );
		if ( ClusterPool::Share()->GetCluster( neighCluster, _para.energyType ) )
		{
			Clusters &one = _pop[m];
			InitTool::AtomTypeSwap( one, neighCluster, _para.energyType );
			m ++;
		} 
	}


	double r0 = AtomPara::NearestNeighborSeparation(_para.atomTypes);
	double length = r0 * pow(_para.N * 1.0, 1.0/3);
	for (int i = m; i < _para.popSize; i++)
	{
		BCGAClusters& cluster = _pop[i];
		InitTool::RandInCubic(cluster,length);
		InitTool::RandType(cluster);
		this->Evaluate(cluster);
	}
}

#pragma endregion 不同的初始化

#pragma region 选择

vector<BCGAClusters> BCGAAlgorithm::Selection()
{
	int popSize = _para.popSize;
	int Noff = (int)(popSize * _para.crossoverRate);
	vector<double> f(popSize);

	this->Fitness();

	double sumOfFitness = 0;
	for ( int i = 0; i < popSize; i++ )
	{
		f[i] = _pop[i].fitness;
		sumOfFitness += f[i];
	}
	for (int i = 0; i < popSize; i++)
		f[i] = f[i] / sumOfFitness;

	for (int i = 1; i < popSize; i++)
		f[i] = f[i-1] + f[i];

	vector<BCGAClusters> seletedClusters;
	int i = 0,lastIndex=-1,index=0;
	while( i < Noff)
	{
		double randInOne = RANDI;
		index = 0;
		while (randInOne > f[index]) index++;
		if ( lastIndex == index ) continue;
		seletedClusters.push_back(_pop[index]);
		lastIndex = index;
		i ++;
	}

	return seletedClusters;
}

void BCGAAlgorithm::Fitness()
{
	double maxEnergy = -DBL_MAX;
	double minEnergy = DBL_MAX;
	int popSize = _para.popSize;

	for (int i =0 ; i < popSize; i ++ )
	{
		double energyOfClusterAtI = _pop[i].GetEnergy();
		if (energyOfClusterAtI > maxEnergy) maxEnergy = energyOfClusterAtI;
		if (energyOfClusterAtI < minEnergy) minEnergy = energyOfClusterAtI;
	}

	for (int i = 0; i < popSize; i++)
	{	
		if ( maxEnergy == minEnergy)
			_pop[i].fitness = 1;
		else
			_pop[i].fitness = (_pop[i].GetEnergy() - minEnergy) / (maxEnergy - minEnergy);
		_pop[i].fitness = exp(-3 * _pop[i].fitness);
	}

	//sort(_pop.begin(),_pop.end(),SortClusterByEnergy);
	//for (int i = 0; i < _pop.size(); i++ )
	//{
	//	cout<<_pop[i].GetEnergy()<<"\t"<<_pop[i].fitness<<endl;
	//}
}	

#pragma endregion 选择

#pragma region 交叉

vector<BCGAClusters> BCGAAlgorithm::CrossOver(vector<BCGAClusters>& parents)
{
	vector<BCGAClusters> offSpringClusters;

	switch (_para.crossover)
	{
	case BCGA_Crossover_Value:
		this->ValueCrossOver( offSpringClusters, parents );
		break;
	case BCGA_Crossover_PlaneCutSplice:
		this->PlaneCutSpliceCrossOver( offSpringClusters, parents );
		break;
	default:
		this->PlaneCutSpliceCrossOver( offSpringClusters, parents );
	}

	return offSpringClusters;
}

void BCGAAlgorithm::ValueCrossOver( vector<BCGAClusters> &children, vector<BCGAClusters> &parents )
{
	int size = parents.size();
	if ( size % 2 == 1 )
		size --;

	for (int i = 0; i < size; i+=2 )
	{
		BCGAClusters fr = parents[i], mr = parents[i+1];
		int cutPoint = RANDIUINT( 3*_para.N-1 );
		
		int atomNumber = (cutPoint+1) / 3;
		for ( int i = 0; i < atomNumber; i++ )
		{
			AtomPos pos1 = fr.GetAtomAtIndex( i ).GetPos();
			AtomPos pos2 = mr.GetAtomAtIndex( i ).GetPos();

			fr.SetPosOfAtomIndex( i, pos2 );
			mr.SetPosOfAtomIndex( i, pos1 );
		}

		int n = atomNumber;
		AtomPos pos1 = fr.GetAtomAtIndex( n ).GetPos();
		AtomPos pos2 = mr.GetAtomAtIndex( n ).GetPos();
		int coordNumber = (cutPoint + 1) % 3;
		if ( coordNumber >= 1 )
		{
			double temp = pos1.x;
			pos1.x = pos2.x;
			pos2.x = temp;
		}
		if ( coordNumber >= 2 )
		{
			double temp = pos1.y;
			pos1.y = pos2.y;
			pos2.y = temp;
		}
		fr.SetPosOfAtomIndex( n, pos1 );
		mr.SetPosOfAtomIndex( n, pos2 );

		this->Evaluate(fr); this->Evaluate(mr);
		children.push_back(fr);
		children.push_back(mr);
	}
}

void BCGAAlgorithm::PlaneCutSpliceCrossOver( vector<BCGAClusters> &children, vector<BCGAClusters> &parents )
{
	int size = parents.size();
	if ( size % 2 == 1 )
		size --;

	for (int i = 0; i < size ; i += 2 )
	{
		BCGAClusters fr = parents[i], mr = parents[i+1];
		CrossOverTool::PlaneCutSplice(fr,mr);
		this->Evaluate(fr); this->Evaluate(mr);
		children.push_back(fr);
		children.push_back(mr);
	}
}

#pragma endregion 交叉

#pragma region 变异

vector<BCGAClusters> BCGAAlgorithm::Mutation(vector<BCGAClusters>& preMutationClusters)
{
	vector<BCGAClusters> mutatedClusters;

	vector<BCGA_Mutate> mutateStrategy;
	if ( (_para.mutate & BCGA_Mutate_Displace) != 0 )
		mutateStrategy.push_back( BCGA_Mutate_Displace );
	if ( (_para.mutate & BCGA_Mutate_Twitst) != 0 )
		mutateStrategy.push_back( BCGA_Mutate_Twitst );
	if ( (_para.mutate & BCGA_Mutate_Replace) != 0 )
		mutateStrategy.push_back( BCGA_Mutate_Replace );
	if ( (_para.mutate & BCGA_Mutate_Swap) != 0 )
		mutateStrategy.push_back( BCGA_Mutate_Swap );

	for (int i = 0; i < preMutationClusters.size(); i++)
	{
		if ( RANDI < _para.muteRate )
		{
			BCGAClusters mutCluster = preMutationClusters[i];
			//cout<<mutCluster.GetEnergy()<<endl;

			int mutateIndex = RANDIUINT( mutateStrategy.size() );
			switch ( mutateStrategy[mutateIndex] )
			{
			case BCGA_Mutate_Displace: MutateTool::AtomDisplacement( mutCluster, _para.mutateNum ); break;
			case BCGA_Mutate_Twitst: MutateTool::Twisting( mutCluster ); break;
			case BCGA_Mutate_Replace: MutateTool::Replacement( mutCluster ); break;
			case BCGA_Mutate_Swap: MutateTool::AtomRandSwap( mutCluster, _para.mutateNum ); break;
			default: MutateTool::Replacement( mutCluster ); break;
			}

			//AtomPos pos = mutCluster.GetAtomAtIndex(0).GetPos();
			//mutCluster.SetPosOfAtomIndex( 1, pos );
			//this->_pe->ForceValue(mutCluster);

			BCGAClusters temp = mutCluster;
			this->Evaluate(mutCluster);
			if ( _finite(mutCluster.GetEnergy()) != 0 )
			{
				//for (int j = 0;j< mutCluster.GetAtomsNumber(); j++ )
				//{
				//	AtomPos pos = temp.GetAtomAtIndex(j).GetPos();
				//	cout<<pos.x<<" "<<pos.y<<" "<<pos.z<<endl;
				//}
				//for (int j = 0;j< mutCluster.GetAtomsNumber(); j++ )
				//{
				//	AtomPos pos = mutCluster.GetAtomAtIndex(j).GetPos();
				//	cout<<pos.x<<" "<<pos.y<<" "<<pos.z<<endl;
				//}
				mutatedClusters.push_back( mutCluster );
			}
			//cout<<mutCluster.GetEnergy()<<endl;

		}
	}

	return mutatedClusters;
}


#pragma endregion 变异

#pragma region 自然选择

void BCGAAlgorithm::Choose(vector<BCGAClusters>& mat, vector<BCGAClusters>& mut)
{
	sort(_pop.begin(),_pop.end(),SortClusterByEnergy);
	sort(mat.begin(),mat.end(),SortClusterByEnergy);
	sort(mut.begin(),mut.end(),SortClusterByEnergy);

	//for ( int i = 1; i < _pop.size(); i++ )
	//{
	//	if ( i < _pop.size() ) cout<<_pop[i].GetEnergy();
	//	cout<<"\t";
	//	if ( i < mat.size() ) cout<<mat[i].GetEnergy();
	//	cout<<"\t";
	//	if ( i < mut.size() ) cout<<mut[i].GetEnergy();
	//	cout<<endl;
	//}
	int popSize = _pop.size();
	int matSize = mat.size();
	int mutSize = mut.size();

	for (int i = matSize - 1; i >= 0; i--)
	{
		BCGAClusters& cluster = mat[i];
		if ( cluster.GetEnergy() > _pop[0].GetEnergy() )
			break;

		int j;
		for (j = 1; j < popSize && _pop[j].GetEnergy() > cluster.GetEnergy(); j++)
			_pop[j-1] = _pop[j];
		_pop[j-1] = cluster;
	}

	for (int i = mutSize-1; i >= 0; i--)
	{
		BCGAClusters& cluster = mut[i];
		if ( cluster.GetEnergy() > _pop[0].GetEnergy() )
			break;

		int j;
		for (j = 1; j < popSize && _pop[j].GetEnergy() > cluster.GetEnergy(); j++)
			_pop[j-1] = _pop[j];
		_pop[j-1] = cluster;
	}

	//for ( int i = 1; i < _pop.size(); i++ )
	//	cout<<_pop[i].GetEnergy()<<endl;

}

void BCGAAlgorithm::Choose1(vector<BCGAClusters>& mat, vector<BCGAClusters>& mut)
{
	sort(_pop.begin(),_pop.end(),SortClusterByEnergy);
	sort(mat.begin(),mat.end(),SortClusterByEnergy);

	int popSize = _pop.size();
	int matSize = mat.size();
	int mutSize = mut.size();

	for (int i = matSize - 1; i >= 0; i--)
	{
		BCGAClusters& cluster = mat[i];
		if ( cluster.GetEnergy() > _pop[0].GetEnergy() )
			break;

		int j;
		for (j = 1; j < popSize && _pop[j].GetEnergy() > cluster.GetEnergy(); j++)
			_pop[j-1] = _pop[j];
		_pop[j-1] = cluster;
	}

	for (int i = 0; i < mutSize; i++)
	{
		_pop[i] = mut[i];
	}


	sort(_pop.begin(),_pop.end(),SortClusterByEnergy);
}

#pragma endregion 自然选择

#pragma region 数据存储

void BCGAAlgorithm::SaveBestCluster()
{
	Clusters &bestCluster = _pop[_pop.size()-1];
	double poolEnergy =	ClusterPool::Share()->GetEnergy( _para.energyType, _para.atomTypes, _para.atomNumber, _para.N);

	if ( bestCluster.GetEnergy() < poolEnergy )
		ClusterPool::Share()->SaveCluster( bestCluster, _para.energyType );
}

void BCGAAlgorithm::SaveEnergyFile()
{
	string folderName = _para.rootFolderName;
	FileTool::CreatFolder( folderName );
	string energyFileName = FileTool::FileNameWithType( folderName, File_Energy );

	double aveEnergy = 0;
	for (int i = 0; i < _pop.size(); i++ )
		aveEnergy += _pop[i].GetEnergy();
	aveEnergy /= _pop.size();

	ofstream file( energyFileName, ios::app );
	file<<setw(10)<<setiosflags(ios::left)<<_iteration;
	file<<setw(20)<<setiosflags(ios::left)<<_pop[_pop.size()-1].GetEnergy();
	file<<setw(20)<<setiosflags(ios::left)<<_pop[0].GetEnergy();
	file<<setw(20)<<setiosflags(ios::left)<<aveEnergy<<endl;
	file.close();
}

#pragma endregion 数据存储

#pragma region 其他

void BCGAAlgorithm::Evaluate( Clusters& cluster )
{
	double E = LocalTool::LocalMinimizeClusters(cluster,_para.energyType);
	cluster.SetEnergy(E);
}

void BCGAAlgorithm::PrintfEnergy( vector<BCGAClusters>& clusters )
{
	cout<<"BEGIN"<<endl;
	for (int i = 0; i < clusters.size(); i++ )
	{
		cout<<i<<clusters[i].GetEnergy()<<endl;
	}
	cout<<"END"<<endl;
}

#pragma endregion 其他


