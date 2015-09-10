#ifndef OPTIMIZER_FUNCTION_WRAPPERS_H
#define OPTIMIZER_FUNCTION_WRAPPERS_H

#include <functional>
#include "OptimizerFunction.hpp"
#include "VAGrisFunction.h"
#include "GKLSFunction.h"

class VAGRisFunctionWrapper : public optimizercore::OptimizerFunction
{
private:
	vagrisfunction::VAGrisFunction mFunction;
public:
	VAGRisFunctionWrapper();
	~VAGRisFunctionWrapper();
	void SetFunctionNumber(int num);
	double GetMinXCoordinate() const;
	double GetMinYCoordinate() const;
	double Calculate(const double *x) override;
	double GetMinValue() const;
};

class GKLSFunctionWrapper : public optimizercore::OptimizerFunction
{
private:
	unsigned mDimention;
	double* mTmpArgument;
	gklsfunction::GKLSFunction mFunction;
public:
	GKLSFunctionWrapper();
	~GKLSFunctionWrapper();
	void SetClassType(gklsfunction::GKLSClass type, unsigned);
	void SetFunctionNumber(int num);
	void SetDimention(unsigned value);
	double Calculate(const double *x) override;
	double GetMinXCoordinate();
	double GetMinYCoordinate();
	double GetMinValue() const;
	void GetMinPoint(double *x);
};

#endif