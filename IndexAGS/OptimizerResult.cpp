#include "OptimizerResult.hpp"

using namespace optimizercore;

OptimizerResult::OptimizerResult()
{
	mIsInitialized = false;
}
OptimizerResult::~OptimizerResult()
{}
OptimizerResult::OptimizerResult(const OptimizerSolution& Solution)
{
	mSolution = Solution;
	mIsInitialized = true;
}
OptimizerSolution OptimizerResult::GetSolution() const
{
	CheckIsInitialized();
	return mSolution;
}
void OptimizerResult::CheckIsInitialized() const
{
	if (mIsInitialized == false)
		throw std::exception("Optimizer Solution is not initialized.");
}