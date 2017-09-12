#include "HolderConstEvaluate.hpp"
#include "OptimizerFunctionWrappers.hpp"
#include "CoreUtils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <omp.h>

using namespace optimizercore;

#define NUM_THREADS 4

void RunEvaluationExpOnGrishaginClass(double * constValues, double eps)
{
	SharedVector leftBound = SharedVector(new double[2]);
	SharedVector rightBound = SharedVector(new double[2]);
	leftBound.get()[0] = leftBound.get()[1] = 0;
	rightBound.get()[0] = rightBound.get()[1] = 1;

	VAGRisFunctionWrapper **functions = utils::AllocateMatrix<VAGRisFunctionWrapper>(NUM_THREADS, 1);
	
	OptimizerMap map(2, 12, MapType::Simple);
	OptimizerSpaceTransformation transform(leftBound, rightBound, 2);
	double holderConsts[100];
	char filename[100];

#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 1; i <= 100; i++)
	{
		functions[omp_get_thread_num()]->SetFunctionNumber(i);
		holderConsts[i - 1] = EvaluateHolderConst(functions[omp_get_thread_num()], map, transform, eps);
		printf("Evaluation for function %i completed. HC = %f\n", i, holderConsts[i - 1]);
	}

	sprintf(filename, "VAGrisFunctions eps=%f.csv", eps);
	FILE* output = fopen(filename, "wt");
	fprintf(output, "Evaluation of holder consts for Grishagin functions at grid with step = %f\n", eps);
	for (int i = 1; i <= 100; i++)
		fprintf(output, "%i; %f\n", i, holderConsts[i - 1]);

	fclose(output);
	utils::DeleteMatrix(functions, NUM_THREADS);
}

void RunEvaluationExpOnGKLSClass(double* constValues, double eps)
{
	int gklsDimension = 2;

	SharedVector leftBound = SharedVector(new double[gklsDimension]);
	SharedVector rightBound = SharedVector(new double[gklsDimension]);
	std::fill_n(leftBound.get(), gklsDimension, -1.);
	std::fill_n(rightBound.get(), gklsDimension, 1.);

	GKLSFunctionWrapper **functions = utils::AllocateMatrix<GKLSFunctionWrapper>(NUM_THREADS, 1);
	for (int i = 0; i < NUM_THREADS; i++)
		functions[i]->SetClassType(gklsfunction::GKLSClass::Simple, gklsDimension);

	double holderConsts[100];
	char filename[100];

	OptimizerMap map(gklsDimension, 12, MapType::Simple);
	OptimizerSpaceTransformation transform(leftBound, rightBound, gklsDimension);

#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 1; i <= 100; i++)
	{
		functions[omp_get_thread_num()]->SetFunctionNumber(i);
		holderConsts[i - 1] = EvaluateHolderConst(functions[omp_get_thread_num()], map, transform, eps);
		printf("Evaluation for function %i completed. HC = %f\n", i, holderConsts[i - 1]);
	}

	sprintf(filename, "GKLSFunctions eps=%f.csv", eps);
	FILE* output = fopen(filename, "wt");
	fprintf(output, "Evaluation of holder consts for GKLS functions at grid with step = %f\n", eps);

	for (int i = 1; i <= 100; i++)
		fprintf(output, "%i; %f\n", i, holderConsts[i - 1]);

	fclose(output);
	utils::DeleteMatrix(functions, NUM_THREADS);
}

double EvaluateHolderConst(OptimizerFunction *function, OptimizerMap& map,
	OptimizerSpaceTransformation& transform, double eps)
{
	int iterationsNumber = (int)(1.0 / eps);
	double holderConst = 0, x = 0;// , normEps = pow(eps, 1.0 / transform.GetDomainDimension());
	double *arg = new double[transform.GetDomainDimension()];
	double *fvalues = new double[iterationsNumber];
	double powValue = 1.0 / transform.GetDomainDimension();

	map.GetImage(x, arg);
	transform.Transform(arg, arg);

	for (int i = 0; i < iterationsNumber - 1; i++) {
		fvalues[i] = function->Calculate(arg);
		for (int j = 0; j < i; j++)
			holderConst = fmax(holderConst, fabs(fvalues[i] - fvalues[j]) /
				pow(eps*(i - j), powValue));
		x += eps;
		map.GetImage(x, arg);
		transform.Transform(arg, arg);
	}

	delete[] arg;
	delete[] fvalues;
	return holderConst;
}
