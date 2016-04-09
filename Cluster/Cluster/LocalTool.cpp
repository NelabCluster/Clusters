#include "StdAfx.h"
#include <lbfgs.h>
#include "LocalTool.h"


static Clusters* localCluster;
static PotentialEnergy* localPE;

static lbfgsfloatval_t evaluate(
	void *instance,
	const lbfgsfloatval_t *x,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
	)
{
	int i;
	lbfgsfloatval_t fx = 0.0;
	int N = n /3;

	for (i = 0; i < N; i++)
		localCluster->SetPosOfAtomIndex(i,AtomPos(x[i],x[N + i],x[2*N + i]));

	fx = localPE->EnergyValue(*localCluster);
	localPE->ForceValue(*localCluster);
	for (i = 0 ;i < N; i++)
	{
		g[i] =  localCluster->GetForceXOfAtomIndex(i);
		g[N + i] =  localCluster->GetForceYOfAtomIndex(i);
		g[2*N + i] =  localCluster->GetForceZOfAtomIndex(i);
	}

	return fx;
}

static int progress(
	void *instance,
	const lbfgsfloatval_t *x,
	const lbfgsfloatval_t *g,
	const lbfgsfloatval_t fx,
	const lbfgsfloatval_t xnorm,
	const lbfgsfloatval_t gnorm,
	const lbfgsfloatval_t step,
	int n,
	int k,
	int ls
	)
{
	/*printf("Iteration %d:\n", k);
	printf("  fx = %f\n", fx);
	printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
	printf("\n");*/
	return 0;
}

LocalTool::LocalTool(void)
{
}


LocalTool::~LocalTool(void)
{
}

double LocalTool::LocalMinimize(Clusters& beforeLocal, Clusters& afterLocal,PE_TYPE type)
{
	afterLocal = beforeLocal;
	return LocalTool::LocalMinimizeClusters(afterLocal,type);
}


double LocalTool::LocalMinimizeClusters(Clusters& cluster,PE_TYPE type)
{
	localCluster = &cluster;
	localPE = PotentialEnergy::initWithType(type);

	int N = cluster.GetAtomsNumber();
	double *cood;
	
	cood = new double[N*3];
	for (int i = 0; i < N; i++)
	{
		AtomPos pos = cluster.GetAtomAtIndex(i).GetPos();
		cood[i] = pos.x; cood[N + i] = pos.y; cood[2*N+i] = pos.z;
	}

	int ret;
	lbfgsfloatval_t fx;
	lbfgs_parameter_t param;

	lbfgs_parameter_init(&param);
	ret = lbfgs(N*3, cood, &fx,evaluate, progress, NULL, &param);


	delete[] cood;
	delete localPE;
	localCluster = NULL;

	return fx;
}



