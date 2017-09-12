#include "OptimizerSpaceTransformation.hpp"
#include "OptimizerSTLFunctionWrapper.hpp"

#include <cmath>

using namespace optimizercore;

#define CUBE_DETECT_PRECISION 0.000000001

OptimizerSpaceTransformation::OptimizerSpaceTransformation() : mIsInitialized(false)
{}

OptimizerSpaceTransformation::OptimizerSpaceTransformation(
	SharedVector leftBound, SharedVector rightBound,
	int domainDimention)
{
	mLeftDomainBound = leftBound;
	mRightDomainBound = rightBound;
	mDimension = domainDimention;
	mLeftBoundPtr = mLeftDomainBound.get();
	mRightBoundPtr = mRightDomainBound.get();

	mRho = 0;
	mSpaceShiftValues = SharedVector(new double[mDimension], utils::array_deleter<double>());
	mSpaceShiftValuesPtr = mSpaceShiftValues.get();

	for (int i = 0; i < mDimension; i++)
	{
		mRho = fmax(mRho, mRightBoundPtr[i] - mLeftBoundPtr[i]);
		mSpaceShiftValuesPtr[i] = 0.5*(mLeftBoundPtr[i] + mRightBoundPtr[i]);
	}

	mNeedZeroConstraint = false;
	for (int i = 0; i < mDimension; i++)
		if (fabs(mRho - (mRightBoundPtr[i] - mLeftBoundPtr[i])) > CUBE_DETECT_PRECISION)
		{
			mNeedZeroConstraint = true;
			break;
		}

	int dimension = mDimension;
	double* spaceShift = mSpaceShiftValuesPtr;
	double* leftBnd = mLeftBoundPtr;
	double* rightBnd = mRightBoundPtr;
	double rho = mRho;

	mZeroConstraint = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
		std::function<double(const double *)>(
		[dimension, spaceShift, rho, leftBnd, rightBnd](const double* x)->double
	{
		double value = -HUGE_VAL;

		for (int i = 0; i < dimension; i++)
			value = fmax(value, (fabs(x[i] - spaceShift[i])
			- (rightBnd[i] - leftBnd[i])*0.5) / rho);

		return value;
	}
	)));

	mIsInitialized = true;
}

void OptimizerSpaceTransformation::Transform(
	const double* arg, double* image) const
{
	for (int i = 0; i < mDimension; i++)
		image[i] = mRho*arg[i] + mSpaceShiftValuesPtr[i];
}

void OptimizerSpaceTransformation::InvertTransform(
	const double* arg, double* image) const
{
	for (int i = 0; i < mDimension; i++)
		image[i] = (arg[i] - mSpaceShiftValuesPtr[i]) / mRho;
}

SharedVector OptimizerSpaceTransformation::GetLeftDomainBound() const
{
	return mLeftDomainBound;
}

SharedVector OptimizerSpaceTransformation::GetRightDomainBound() const
{
	return mRightDomainBound;
}

int OptimizerSpaceTransformation::GetDomainDimension() const
{
	return mDimension;
}
OptimizerFunctionPtr OptimizerSpaceTransformation::GetZeroConstraint() const
{
	return mZeroConstraint;
}
bool OptimizerSpaceTransformation::IsZeroConstraintActive() const
{
	return mNeedZeroConstraint;
}