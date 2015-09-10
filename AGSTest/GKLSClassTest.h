#ifndef GKLS_CLASS_TEST_H
#define GKLS_CLASS_TEST_H
#include "OptimizerFunctionWrappers.h"
#include "OptimizerAlgorithm.hpp"

int ParseArguments(int arg_c, char** arg_v, int &threadsNum, int &problemDim,
	gklsfunction::GKLSClass& classType, int& maxIterCount);
void TestGKLSClass(optimizercore::OptimizerParameters, gklsfunction::GKLSClass, int gklsDimention);

#endif