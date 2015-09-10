#include "OptimizerSolution.hpp"
#include <cassert>

using namespace optimizercore;

OptimizerSolution::OptimizerSolution()
{
	mIsInitialized = false;
}
OptimizerSolution::OptimizerSolution(int iterationsCount, double optimumValue,
	double oneDimOptimumPoint, unsigned dimention, SharedVector optimumPoint)
{
	assert(iterationsCount > 0);
	assert(dimention > 0);

	mIteratinsCount = iterationsCount;
	mOptimumValue = optimumValue;
	mOneDimOptimumPoint = oneDimOptimumPoint;
	mOptimumPoint = optimumPoint;
	mDimention = dimention;

	mIsInitialized = true;
}
OptimizerSolution::~OptimizerSolution()
{

}

int OptimizerSolution::GetIterationsCount() const
{
	CheckIsInitialized();
	return mIteratinsCount;
}
double OptimizerSolution::GetOneDimOptimumPoint() const
{
	CheckIsInitialized();
	return mOneDimOptimumPoint;
}
double OptimizerSolution::GetOptimumValue() const
{
	CheckIsInitialized();
	return mOptimumValue;
}
SharedVector OptimizerSolution::GetOptimumPoint() const
{
	CheckIsInitialized();
	return mOptimumPoint;
}

void OptimizerSolution::CheckIsInitialized() const
{
	if (mIsInitialized == false)
		throw std::exception("Optimizer Solution is not initialized.");
}