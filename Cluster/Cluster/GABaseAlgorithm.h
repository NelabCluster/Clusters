#pragma once
#include "ialgorithm.h"
#include "GABaseParameters.h"
#include "GABaseClusters.h"

class GABaseAlgorithm :
	public IAlgorithm
{
public:
	GABaseAlgorithm(GABaseParameters &para);
	~GABaseAlgorithm(void);

	void Initialization(void);
	void Process(void);
	bool EndCondition(void);

private:
	GABaseParameters _parameters;
	//GABaseClusters _bestClusters;
	vector<GABaseClusters> _pop;

	vector<GABaseClusters> Selection();
	void CrossOver();
	void Mutation();
};

