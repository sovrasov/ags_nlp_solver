#include "OptimizerSTLFunctionWrapper.hpp"

using namespace optimizercore;

OptimizerSTLFunctionWrapper::OptimizerSTLFunctionWrapper(std::function<double(const double*)> function)
{
	mFunction = function;
	mCalculationsCounter = 0;
}
void OptimizerSTLFunctionWrapper::SetFunction(std::function<double(const double*)> function)
{
	mFunction = function;
	mCalculationsCounter = 0;
}
double OptimizerSTLFunctionWrapper::Calculate(const double *x)
{
	mCalculationsCounter++;
	return mFunction(x);
}
int OptimizerSTLFunctionWrapper::GetCalculationsCounter() const
{
	return mCalculationsCounter;
}
void OptimizerSTLFunctionWrapper::ResetCalculationsCounter()
{
	mCalculationsCounter = 0;
}