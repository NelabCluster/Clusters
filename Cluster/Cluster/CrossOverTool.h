#pragma once
#include "stdafx.h"
#include "Clusters.h"

class CrossOverTool
{
public:
	static void PlaneCutSpliceOperate(Clusters& top,Clusters& bottom,int cutPoint);
	static void PlaneCutSplice(Clusters& mother, Clusters& father);
	
};

