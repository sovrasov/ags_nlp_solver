#include "OptimizerMap.hpp"
#include <cassert>

optimizercore::OptimizerMap::OptimizerMap()
{
}

optimizercore::OptimizerMap::~OptimizerMap()
{
}

optimizercore::OptimizerMap::OptimizerMap(int dimension, int tightness, MapType type)
{
	assert(dimension > 1);
	assert(tightness > 2);

	mDimension = dimension;
	mTightness = tightness;
	mMapType = type;
}

void optimizercore::OptimizerMap::GetImage(double x, double y[])
{
	switch (mMapType)
	{
	case MapType::Simple:
		mapd(x, mTightness, y, mDimension, 1);
		break;
	case MapType::Linear:
		mapd(x, mTightness, y, mDimension, 2);
		break;
	case MapType::Noninjective:
		mapd(x, mTightness, y, mDimension, 3);
		break;
	default:
		mapd(x, mTightness, y, mDimension, 1);
	}
}

int optimizercore::OptimizerMap::GetAllPreimages(double * p, double xp[])
{
	int preimNumber = 1;
	switch (mMapType)
	{
	case MapType::Simple:
		xyd(xp, mTightness, p, mDimension);
		break;
	case MapType::Linear:
		xyd(xp, mTightness, p, mDimension);
		break;
	case MapType::Noninjective:
		invmad(mTightness, xp, MAX_PREIMAGES, &preimNumber, p, mDimension, 4);
		break;
	default:
		xyd(xp, mTightness, p, mDimension);
	}
	return preimNumber;
}
