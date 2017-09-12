#include "OptimizerMap.hpp"
#include <cassert>

optimizercore::OptimizerMap::OptimizerMap()
{
	mIsInitialized = false;
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

	switch (mMapType)
	{
	case MapType::Simple:
		mMapKey = 1;
		break;
	case MapType::Linear:
		mMapKey = 2;
		break;
	case MapType::Noninjective:
		mMapKey = 3;
		break;
	}

	mIsInitialized = true;
}

void optimizercore::OptimizerMap::GetImage(double x, double y[])
{
	mapd(x, mTightness, y, mDimension, mMapKey);
}

int optimizercore::OptimizerMap::GetAllPreimages(double * p, double xp[])
{
	int preimNumber = 1;
	if(mMapType == MapType::Noninjective)
		invmad(mTightness, xp, MAX_PREIMAGES, &preimNumber, p, mDimension, 4);
	else
		xyd(xp, mTightness, p, mDimension);

	return preimNumber;
}
