#include "OptimizerTask.hpp"
#include <cassert>

using namespace optimizercore;

OptimizerTask::~OptimizerTask()
{}
OptimizerTask::OptimizerTask() : mIsInitialized(false)
{}
OptimizerTask::OptimizerTask(std::shared_ptr<OptimizerFunctionPtr> taskFunctions,
	unsigned numberOfRestrictions, unsigned taskDimention,
	SharedVector leftBound, SharedVector rightBound)
{
	assert(taskDimention > 0);
	assert(numberOfRestrictions >= 0);
	assert(taskFunctions);
	assert(leftBound);
	assert(rightBound);

	mTaskFunctions = taskFunctions;
	mTaskDimension = taskDimention;
	mNumberOfRestrictions = numberOfRestrictions;

	mSpaceTransform = OptimizerSpaceTransformation(leftBound, rightBound, mTaskDimension);

	mOptimumPoint = nullptr;
	mIsInitialized = true;
}
unsigned OptimizerTask::GetNumberOfRestrictions() const
{
	CheckIsInitialized();
	return mNumberOfRestrictions;
}
unsigned OptimizerTask::GetTaskDimension() const
{
	CheckIsInitialized();
	return mTaskDimension;
}
SharedVector OptimizerTask::GetOptimumPoint() const
{
	CheckIsInitialized();
	return mOptimumPoint;
}
void OptimizerTask::SetOptimumPoint(SharedVector value)
{
	mOptimumPoint = value;
}
std::shared_ptr<OptimizerFunctionPtr> OptimizerTask::GetTaskFunctions() const
{
	CheckIsInitialized();
	return mTaskFunctions;
}
OptimizerSpaceTransformation OptimizerTask::GetSpaceTransformation() const
{
	return mSpaceTransform;
}
void OptimizerTask::CheckIsInitialized() const
{
	if (mIsInitialized == false)
		throw std::runtime_error("Optimizer Task is not initialized.");
}