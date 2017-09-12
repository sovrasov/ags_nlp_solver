#include "OptimizerResult.hpp"

using namespace optimizercore;

OptimizerResult::OptimizerResult()
{
	mIsInitialized = false;
}
optimizercore::OptimizerResult::OptimizerResult(const OptimizerSolution & Solution)
{
	mSolution = Solution;
	mNumberOfFunctionals = 1;
	mFunctionalsCalculationStat = 0;
	mIsInitialized = true;
}
OptimizerResult::~OptimizerResult()
{}
OptimizerResult::OptimizerResult(const OptimizerSolution& Solution,
	SharedIntVector functionalsCalculationStat,
	int numberOfFunctionals)
{
	mSolution = Solution;
	mFunctionalsCalculationStat = functionalsCalculationStat;
	mNumberOfFunctionals = numberOfFunctionals;

	mIsInitialized = true;
}
OptimizerSolution OptimizerResult::GetSolution() const
{
	CheckIsInitialized();
	return mSolution;
}
int optimizercore::OptimizerResult::GetNumberOfCalculations(int fNumber) const
{
	return mFunctionalsCalculationStat.get()[fNumber];
}
int optimizercore::OptimizerResult::GetNumberOfFunctionals() const
{
	return mNumberOfFunctionals;
}
void OptimizerResult::CheckIsInitialized() const
{
	if (mIsInitialized == false)
		throw std::runtime_error("Optimizer Solution is not initialized.");
}
