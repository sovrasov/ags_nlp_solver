#include "HookeJeevesLocalMethod.hpp"
#include <cassert>
#include <cstring>
#include <cmath>
#include <algorithm>

#define MAX_LOCAL_ITERATIONS_NUMBER 10000

using namespace localoptimizer;

HookeJeevesLocalMethod::HookeJeevesLocalMethod()
{
	mCurrentPoint = nullptr;
	mPreviousResearchDirection = nullptr;
	mCurrentResearchDirection = nullptr;
	mStartPoint = nullptr;
	mFunctions = nullptr;

	mDimension = 0;

	mEps = 0.001;
	mStep = 0.01;
	mStepMultiplier = 2;
}
HookeJeevesLocalMethod::~HookeJeevesLocalMethod()
{
	if (mFunctions)
		delete[] mFunctions;
	if (mStartPoint)
		delete[] mStartPoint;
}
void HookeJeevesLocalMethod::DoStep()
{
	for (int i = 0; i < mDimension; i++)
		mCurrentPoint[i] = (1 + mStepMultiplier)*mCurrentResearchDirection[i] -
			mStepMultiplier*mPreviousResearchDirection[i];
}
void HookeJeevesLocalMethod::StartOptimization(double* optimumPointEvaluation)
{
	int k = 0, i=0;
	bool needRestart = true;
	double currentFValue, nextFValue;

	mCurrentPoint = new double[mDimension];
	mCurrentResearchDirection = new double[mDimension];
	mPreviousResearchDirection = new double[mDimension];

	while (i < MAX_LOCAL_ITERATIONS_NUMBER)	{
		i++;
		if (needRestart)	{
			k = 0;
			std::memcpy(mCurrentPoint, mStartPoint, sizeof(double)*mDimension);
			std::memcpy(mCurrentResearchDirection, mStartPoint, sizeof(double)*mDimension);
			currentFValue = EvaluateTargetFunctiuon(mCurrentPoint);
			needRestart = false;
		}

		std::swap(mPreviousResearchDirection, mCurrentResearchDirection);
		std::memcpy(mCurrentResearchDirection, mCurrentPoint, sizeof(double)*mDimension);
		nextFValue = MakeResearch(mCurrentResearchDirection);

		if (currentFValue > nextFValue)	{
			DoStep();
			k++;
			currentFValue = nextFValue;
		}
		else if (mStep > mEps)	{
			if (k != 0)
				std::swap(mStartPoint, mPreviousResearchDirection);
			else
				mStep /= mStepMultiplier;
			needRestart = true;
		}
		else
			break;
	}

	std::memcpy(optimumPointEvaluation, mPreviousResearchDirection, sizeof(double)*mDimension);

	delete[] mCurrentPoint;
	delete[] mPreviousResearchDirection;
	delete[] mCurrentResearchDirection;
}
double HookeJeevesLocalMethod::EvaluateTargetFunctiuon(const double* x) const
{
	for (int i = 0; i <= mConstraintsNumber; i++)
	{
		double value = mFunctions[i]->Calculate(x);
		if (i < mConstraintsNumber && value > 0)
			return HUGE_VAL;
		else if (i == mConstraintsNumber)
			return value;
	}
	return HUGE_VAL;
}
void HookeJeevesLocalMethod::SetEps(double eps)
{
	assert(eps > 0);
	mEps = eps;
}
void HookeJeevesLocalMethod::SetInitialStep(double value)
{
	assert(value > 0);
	mStep = value;
}
void HookeJeevesLocalMethod::SetStepMultiplier(double value)
{
	assert(mStepMultiplier > 1);
	mStepMultiplier = value;
}
void HookeJeevesLocalMethod::SetProblem(OptimizerTask task)
{
	mConstraintsNumber = task.GetNumberOfRestrictions() + 1;
	if (mFunctions)
		delete[] mFunctions;
	mFunctions = new OptimizerFunction*[mConstraintsNumber + 1];
	mFunctions[0] = task.GetSpaceTransformation().GetZeroConstraint().get();
	for (int i = 1; i <= mConstraintsNumber; i++)
		mFunctions[i] = task.GetTaskFunctions().get()[i - 1].get();
}
void HookeJeevesLocalMethod::SetStartPoint(const double* point, int dimention)
{
	assert(point);
	assert(dimention > 1);
	mDimension = dimention;
	mStartPoint = new double[mDimension];
	std::memcpy(mStartPoint, point, mDimension*sizeof(double));
}
double HookeJeevesLocalMethod::MakeResearch(double* startPoint)
{
	double bestValue = EvaluateTargetFunctiuon(startPoint);

	for (int i = 0; i < mDimension; i++)
	{
		startPoint[i] += mStep;
		double rightFvalue = EvaluateTargetFunctiuon(startPoint);

		if (rightFvalue > bestValue)
		{
			startPoint[i] -= 2 * mStep;
			double leftFValue = EvaluateTargetFunctiuon(startPoint);
			if (leftFValue > bestValue)
				startPoint[i] += mStep;
			else
				bestValue = leftFValue;
		}
		else
			bestValue = rightFvalue;
	}

	return bestValue;
}
