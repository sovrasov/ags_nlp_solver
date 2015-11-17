#include "GKLSClassTest.hpp"
#include "CoreUtils.hpp"
#include "OptimizerAlgorithmUnconstrained.hpp"
#include "OptimizerSpaceTransformation.hpp"

#include <cstdio>
#include <algorithm>
#include <chrono>

using namespace std::chrono;
using namespace optimizercore;

int ParseArguments(int arg_c, char** arg_v, int &threadsNum, int &problemDim,
	gklsfunction::GKLSClass& classType, int& maxIterCount)
{
	if (arg_c > 2)
	{
		for (int i = 1; i < arg_c; i++)	{
			if (!strcmp(arg_v[i], "-numThreads"))
				threadsNum = atoi(arg_v[i + 1]);
			else if (!strcmp(arg_v[i], "-dim"))
				problemDim = atoi(arg_v[i + 1]);
			if (!strcmp(arg_v[i], "-simple"))
				classType = gklsfunction::GKLSClass::Simple;
			if (!strcmp(arg_v[i], "-hard"))
				classType = gklsfunction::GKLSClass::Hard;
			if (!strcmp(arg_v[i], "-iterationsNum"))
				maxIterCount = atoi(arg_v[i + 1]);
		}
		return 0;
	}
	else
		return 1;
}

void TestGKLSClass(optimizercore::OptimizerParameters algParams, gklsfunction::GKLSClass classType, int gklsDimention)
{
	SharedVector leftBound = SharedVector(new double[2]);
	SharedVector rightBound = SharedVector(new double[2]);
	leftBound.get()[0] = leftBound.get()[1] = -1;
	rightBound.get()[0] = rightBound.get()[1] = 1;

	GKLSFunctionWrapper *function = new GKLSFunctionWrapper();
	auto taskFunctions = new optimizercore::OptimizerFunctionPtr[1];
	taskFunctions[0] = OptimizerFunctionPtr(function);

	function->SetClassType(classType, gklsDimention);

	double meanItCount = 0;
	int err_count = 0;
	int results[100];
	std::fill_n(results, 100, 0);
	int max_count=0;
	
//	algParams.eps /= 2;
	//optimizercore::OptimizerAlgorithm ags;
	optimizercore::OptimizerAlgorithmUnconstrained ags;
	ags.SetParameters(algParams);
	optimizercore::OptimizerTask task(std::shared_ptr<OptimizerFunctionPtr>(taskFunctions,
		utils::array_deleter<OptimizerFunctionPtr>()), 0, gklsDimention, leftBound, rightBound);
	//ags.SetTask(task);
	ags.SetTask(taskFunctions[0], 
		optimizercore::OptimizerSpaceTransformation(leftBound, rightBound, gklsDimention));

	double* globalMinPoint = new double[gklsDimention], *y;

	for (int j = 2; j <= 2; j++)
	{
		function->SetFunctionNumber(j);
		function->GetMinPoint(globalMinPoint);
		printf("\nTask# %i\n", j);
		for (unsigned i = 0; i < gklsDimention; i++)
		{
	//		printf("%f  ", globalMinPoint[i]);
		//	globalMinPoint[i] *= 0.5;
		}
	//	printf("\n");
			auto result = ags.StartOptimization(globalMinPoint, optimizercore::StopCriterionType::OptimalPoint);
			auto stat = result.GetSolution();
			y = stat.GetOptimumPoint().get();

			for (unsigned i = 0; i < gklsDimention; i++)
				printf("%f  ", y[i]);

			printf("\nIt_count: %i\n", stat.GetIterationsCount());
			printf("\nFunction value %f\n", stat.GetOptimumValue());
			printf("Helder const evaluation: %f", ags.GetLipschitzConst());
			printf("\n-------------------\n");
			meanItCount += stat.GetIterationsCount() / 100.0;
			if (optimizercore::utils::NormNDimMax(y, globalMinPoint, gklsDimention) > 0.01)
			{
				err_count++;
				results[j - 1] = 21000000;
				//break;
			}
			else{
				results[j - 1] = stat.GetIterationsCount();
				if (max_count < results[j - 1])
					max_count = results[j - 1];
			}
	}

	printf("Total errors: %i\n", err_count);
	printf("Mean iterations number: %f\n", meanItCount);

	if (err_count == 0)
	{
		FILE* out;
		char filename[100];
		int m_type = static_cast<int>(algParams.mapType);

		sprintf_s(filename, "R= %f map_type= %i local percent= %i, threadsNum= %i .txt",
			*algParams.r, m_type, algParams.localMixParameter, algParams.numberOfThreads);



		printf("\nR= %f map_type= %i local percent= %i, threadsNum= %i .txt",
			*algParams.r, m_type, algParams.localMixParameter, algParams.numberOfThreads);

		out = fopen(filename, "wt");
		int t_count = 0;

		for (int i = 0; i < max_count + 20; i += 10){
			for (int j = 0; j < 100; j++)
				if (results[j] < i)
					t_count++;
			fprintf(out, "%i; %i\n", i, t_count);
			//printf("%i; %i\n", i, t_count);
			t_count = 0;
		}
		fclose(out);
	}
}
int TestMultimapsGKLSClass(optimizercore::OptimizerParameters algParams, gklsfunction::GKLSClass classType, int gklsDimention)
{
	algParams.algDimention = gklsDimention;
	SharedVector leftBound = SharedVector(new double[gklsDimention]);
	SharedVector rightBound = SharedVector(new double[gklsDimention]);
	std::fill_n(leftBound.get(), gklsDimention, -1.);
	std::fill_n(rightBound.get(), gklsDimention, 1.);

	GKLSFunctionWrapper *function = new GKLSFunctionWrapper();
	auto taskFunctions = new optimizercore::OptimizerFunctionPtr[1];
	taskFunctions[0] = OptimizerFunctionPtr(function);

	function->SetClassType(classType, gklsDimention);

	double meanItCount = 0;
	int err_count = 0;
	int results[100];
	std::fill_n(results, 100, 0);
	int max_count = 0;

	optimizercore::OptimizerAlgorithm ags;
	ags.SetParameters(algParams);
	optimizercore::OptimizerTask task(std::shared_ptr<OptimizerFunctionPtr>(taskFunctions,
		utils::array_deleter<OptimizerFunctionPtr>()), 0, gklsDimention, leftBound, rightBound);
	ags.SetTask(task);

	double* globalMinPoint = new double[gklsDimention], *y;
	double avgHelderConst = 0;

	steady_clock::time_point start = steady_clock::now();
	for (int j = 1; j <= 100; j++)
	{
		function->SetFunctionNumber(j);
		function->GetMinPoint(globalMinPoint);
		printf("\nTask# %i\n", j);
		for (unsigned i = 0; i < gklsDimention; i++)
		{
			//		printf("%f  ", globalMinPoint[i]);
			//	globalMinPoint[i] *= 0.5;
		}
		//	printf("\n");
		auto result = ags.StartOptimization(globalMinPoint, optimizercore::StopCriterionType::OptimalPoint);
		auto stat = result.GetSolution();
		y = stat.GetOptimumPoint().get();
		double helderConst = ags.GetLipschitzConst(result.GetNumberOfFunctionals() - 1);
		avgHelderConst += helderConst / 100.0;

		for (unsigned i = 0; i < gklsDimention; i++)
			printf("%f  ", y[i]);
		printf("\nIt_count: %i\n", result.GetNumberOfCalculations(0));
		printf("\nFunction value %f\n", stat.GetOptimumValue());
		printf("Helder const evaluation: %f", helderConst);
		printf("\n-------------------\n");
		meanItCount += result.GetNumberOfCalculations(result.GetNumberOfFunctionals() - 1) / 100.0;

		if (optimizercore::utils::NormNDimMax(y, globalMinPoint, gklsDimention) > 0.01)
		{
			err_count++;
			results[j - 1] = 21000000;
			//break;
		}
		else {
			//results[j - 1] = stat.GetIterationsCount();
			results[j - 1] = result.GetNumberOfCalculations(result.GetNumberOfFunctionals() - 1);
			if (max_count < results[j - 1])
				max_count = results[j - 1];
		}
	}
	steady_clock::time_point end = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(end - start);

	printf("Total errors: %i\n", err_count);
	printf("Mean iterations number: %f\n", meanItCount);
	printf("Average helder const: %f\n", avgHelderConst);
	printf("Total time: %f\n", time_span.count());

	if (err_count == 0)
	{
		FILE* out;
		char filename[100];
		int m_type = static_cast<int>(algParams.mapType);

		sprintf_s(filename, "GK R= %.3f map_type= %i local percent= %i, threadsNum= %i maps_num=%i.csv",
			*algParams.r, m_type, algParams.localMixParameter, algParams.numberOfThreads,
			algParams.numberOfMaps);

		printf("\nGK R= %f map_type= %i local percent= %i, threadsNum= %i\n",
			*algParams.r, m_type, algParams.localMixParameter, algParams.numberOfThreads);

		out = fopen(filename, "wt");
		int t_count = 0;

		for (int i = 0; i < max_count + 20; i += 10) {
			for (int j = 0; j < 100; j++)
				if (results[j] < i)
					t_count++;
				fprintf(out, "%i; %i\n", i, t_count);
			//printf("%i; %i\n", i, t_count);
			t_count = 0;
		}
		fprintf(out, "Average helder const: %f\n", avgHelderConst);
		fprintf(out, "Average iterations number: %f\n", meanItCount);
		fclose(out);
	}
	return err_count;
}