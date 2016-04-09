#include "StdAfx.h"
#include "PropertyTool.h"


PropertyTool::PropertyTool(void)
{
}


PropertyTool::~PropertyTool(void)
{
}

void PropertyTool::FingerprintFunction( Clusters &cluster )
{
	int N = cluster.GetAtomsNumber();
	cluster.MoveCenterToOrigin();

	double maxR = 0;
	vector<double> R;
	for ( int i = 0; i < N; i++ )
	{
		AtomPos pos = cluster.GetAtomAtIndex( i ).GetPos();
		double r = sqrt( pos.x * pos.x + pos.y * pos.y + pos.z * pos.z );
		if ( r > maxR)
			maxR = r;
		R.push_back( r );
	}

	double r0 = AtomPara::NearestNeighborSeparation( cluster.GetAlloy() );
	double radius = maxR;//r0 * pow( N * 1.0, 1.0/3 );
	double V = 4 * M_PI * radius * radius * radius;

	vector<int> num;
	vector<double> g;
	double r = 0,dr = 0.1;
	int sum = 0;
	while ( sum != N )
	{
		int temp = 0;
		for ( int i = 0; i < N; i++ )
		{
			if ( R[i] >= r && R[i] < (r+dr)  )
			{
				temp ++;
				sum ++;
			}
		}
		num.push_back( temp );
		double ring = 4 * M_PI * r * r * dr;
		if ( r == 0 )
			ring = 4 * M_PI * dr;
		double gr = V * temp / ring / N;
		g.push_back(gr);
		r += dr;
	}

}