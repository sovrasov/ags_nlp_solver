#ifndef HOLDER_CONST_EVALUATE_HPP
#define HOLDER_CONST_EVALUATE_HPP

#include "OptimizerFunction.hpp"
#include "OptimizerMap.hpp"
#include "OptimizerSpaceTransformation.hpp"

using namespace optimizercore;

void RunEvaluationExpOnGrishaginClass(double* constValues, double eps);
void RunEvaluationExpOnGKLSClass(double* constValues, double eps);

double EvaluateHolderConst(OptimizerFunction *function, OptimizerMap& map,
	OptimizerSpaceTransformation& transform, double eps);

#endif