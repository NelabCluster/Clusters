#pragma once
#include "StdAfx.h"
#include "ialgorithm.h"
#include "BCGAParameters.h"
#include "BCGAClusters.h"

class BCGAAlgorithm :
	public IAlgorithm
{
public:
	BCGAAlgorithm(BCGAParameters& para);
	~BCGAAlgorithm(void);

	void Initialization();
	void Process();
	bool EndCondition();

	double GetBestEnergy();
	void PrintProgress();
	void SaveProgress();
	void Ending();

protected:
	BCGAParameters _para;
	PotentialEnergy* _pe;

	int _iteration;
	int _convergenceCount;
	double _lastPopBestEnergy;

	vector<BCGAClusters> _pop;

	clock_t _beginTime;
	clock_t _endTime;

	void RandInitialization();
	void GrowInitialization();
	void SwapInitialization();

	vector<BCGAClusters> Selection();
	void Fitness();

	vector<BCGAClusters> CrossOver(vector<BCGAClusters>& parents);
	void ValueCrossOver( vector<BCGAClusters> &children, vector<BCGAClusters> &parents );
	void PlaneCutSpliceCrossOver( vector<BCGAClusters> &children, vector<BCGAClusters> &parents );

	vector<BCGAClusters> Mutation(vector<BCGAClusters>& preMutationClusters);

	void Choose(vector<BCGAClusters>& mat, vector<BCGAClusters>& mut);
	void Choose1(vector<BCGAClusters>& mat, vector<BCGAClusters>& mut);

	void SaveBestCluster();
	void SaveEnergyFile();

	void Evaluate(Clusters& cluster);
	void PrintfEnergy( vector<BCGAClusters>& clusters );

};

