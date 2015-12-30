#include "OptimizerExperimentalMap.hpp"

optimizercore::OptimizerExperimentalMap::OptimizerExperimentalMap()
{
}

optimizercore::OptimizerExperimentalMap::~OptimizerExperimentalMap()
{
}

void optimizercore::OptimizerExperimentalMap::GetImage(double x, double y[])
{
	//for now only for 4d points
	double z[2];
	mapd(x, 14, z, 2);
	mapd(z[0] + 0.5, 14, y, 2);
	mapd(z[1] + 0.5, 14, y + 2, 2);
}

int optimizercore::OptimizerExperimentalMap::GetAllPreimages(double * p, double xp[])
{
	return 0;
}
