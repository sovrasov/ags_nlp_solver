#include "OptimizerFunctionWrappers.h"

VAGRisFunctionWrapper::VAGRisFunctionWrapper()
{	}
VAGRisFunctionWrapper::~VAGRisFunctionWrapper()
{	}
void VAGRisFunctionWrapper::SetFunctionNumber(int num)
{
	mFunction.SetFunctionNumber(num);
}
double VAGRisFunctionWrapper::GetMinXCoordinate() const
{
	return mFunction.GetMinimumXCoordinate(mFunction.GetFunctionNumber());
}
double VAGRisFunctionWrapper::GetMinYCoordinate() const
{
	return mFunction.GetMinimumYCoordinate(mFunction.GetFunctionNumber());
}
double VAGRisFunctionWrapper::Calculate(const double *x)
{
	return mFunction.EvaluateFunction(x[0], x[1]);
}
double VAGRisFunctionWrapper::GetMinValue() const
{
	return mFunction.EvaluateFunction(GetMinXCoordinate(), GetMinYCoordinate());
}





GKLSFunctionWrapper::GKLSFunctionWrapper() : mTmpArgument(nullptr), mDimention(2)
{	}
void GKLSFunctionWrapper::SetClassType(gklsfunction::GKLSClass type, unsigned dimention)
{
	mFunction.SetFunctionClass(type, dimention);
	mDimention = mFunction.GetDimention();

	if (mTmpArgument != nullptr)
		delete[] mTmpArgument;
	mTmpArgument = new double[mDimention];
}
void GKLSFunctionWrapper::SetFunctionNumber(int num)
{
	mFunction.SetFunctionNumber(num);
}
double GKLSFunctionWrapper::Calculate(const double *x)
{
	return mFunction.EvaluateDFunction(x);;
}
GKLSFunctionWrapper::~GKLSFunctionWrapper()
{
	if (mTmpArgument != nullptr)
		delete[] mTmpArgument;
}
double GKLSFunctionWrapper::GetMinXCoordinate()
{
	mFunction.GetGlobalMinimumPoint(mFunction.GetFunctionNumber(), mTmpArgument);
	return mTmpArgument[0];
}
double GKLSFunctionWrapper::GetMinYCoordinate()
{
	mFunction.GetGlobalMinimumPoint(mFunction.GetFunctionNumber(), mTmpArgument);
	return mTmpArgument[1];
}
double GKLSFunctionWrapper::GetMinValue() const
{
	return mFunction.GetGlobalMinimumValue();
}
void GKLSFunctionWrapper::GetMinPoint(double *x)
{
	mFunction.GetGlobalMinimumPoint(mFunction.GetFunctionNumber(), x);
}
void GKLSFunctionWrapper::SetDimention(unsigned value)
{
	mDimention = value;
	mFunction.SetDimention(mDimention);
}
