#include "local_optimizer.hpp"

#include <cmath>
#include <algorithm>
#include <limits>

#define MAX_LOCAL_ITERATIONS_NUMBER 10000

void HookeJeevesOptimizer::SetParameters(double eps, double step, double stepMult)
{
  NLP_SOLVER_ASSERT(eps > 0 && step > 0 && stepMult > 0, "Wrong papameters of the local optimizer");
  mEps = eps;
  mStep = step;
  mStepMultiplier = stepMult;
}

Trial HookeJeevesOptimizer::Optimize(std::shared_ptr<IGOProblem<double>> problem, const Trial& startPoint)
{
  mProblem = problem;
  mStartPoint = startPoint;

  int k = 0, i=0;
	bool needRestart = true;
	double currentFValue, nextFValue;

	while (i < MAX_LOCAL_ITERATIONS_NUMBER)	{
		i++;
		if (needRestart)	{
			k = 0;
      mCurrentPoint = mStartPoint;
      mCurrentResearchDirection = mStartPoint;
			//std::memcpy(mCurrentPoint, mStartPoint, sizeof(double)*mDimension);
			//std::memcpy(mCurrentResearchDirection, mStartPoint, sizeof(double)*mDimension);
			currentFValue = ComputeObjective(mCurrentPoint.y);
			needRestart = false;
		}

		std::swap(mPreviousResearchDirection, mCurrentResearchDirection);
    mCurrentResearchDirection = mCurrentPoint;
		//std::memcpy(mCurrentResearchDirection, mCurrentPoint, sizeof(double)*mDimension);
		nextFValue = MakeResearch(mCurrentResearchDirection.y);

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

	//std::memcpy(optimumPointEvaluation, mPreviousResearchDirection, sizeof(double)*mDimension);

  return mPreviousResearchDirection;
}

void HookeJeevesOptimizer::DoStep()
{
	for (int i = 0; i < mProblem->GetDimension(); i++)
		mCurrentPoint.y[i] = (1 + mStepMultiplier)*mCurrentResearchDirection.y[i] -
			mStepMultiplier*mPreviousResearchDirection.y[i];
}

double HookeJeevesOptimizer::ComputeObjective(const double* x) const
{
	for (int i = 0; i <= mProblem->GetConstraintsNumber(); i++)
	{
		double value = mProblem->Calculate(x, i);
		if (i < mProblem->GetConstraintsNumber() && value > 0)
			return std::numeric_limits<double>::max();
		else if (i == mProblem->GetConstraintsNumber())
			return value;
	}
	return std::numeric_limits<double>::max();
}

double HookeJeevesOptimizer::MakeResearch(double* startPoint)
{
	double bestValue = ComputeObjective(startPoint);

	for (int i = 0; i < mProblem->GetDimension(); i++)
	{
		startPoint[i] += mStep;
		double rightFvalue = ComputeObjective(startPoint);

		if (rightFvalue > bestValue)
		{
			startPoint[i] -= 2 * mStep;
			double leftFValue = ComputeObjective(startPoint);
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
/*
	mEps = 0.001;
	mStep = 0.01;
	mStepMultiplier = 2;
*/
