#include "TestGrishaginClass.hpp"
#include "OptimizerCoreGlobal.hpp"

#include "OptimizerFunctionWrappers.hpp"
#include "OptimizerAlgorithm.hpp"

#include <chrono>

using namespace std::chrono;
using namespace optimizercore;

int TestVAGrisClass(optimizercore::OptimizerParameters algParams)
{
	SharedVector leftBound = SharedVector(new double[2]);
	SharedVector rightBound = SharedVector(new double[2]);
	leftBound.get()[0] = leftBound.get()[1] = 0;
	rightBound.get()[0] = rightBound.get()[1] = 1;

	VAGRisFunctionWrapper *function = new VAGRisFunctionWrapper();
	auto taskFunctions = new optimizercore::OptimizerFunctionPtr[1];
	taskFunctions[0] = OptimizerFunctionPtr(function);

	double meanItCount = 0, avgHelderConst = 0;;
	double globalMinPoint[2], *y;
	int err_count = 0;
	int results[100];
	std::fill_n(results, 100, 0);
	int max_count = 0;

	optimizercore::OptimizerAlgorithm ags;
	ags.SetParameters(algParams);
	optimizercore::OptimizerTask task(std::shared_ptr<OptimizerFunctionPtr>(taskFunctions,
		utils::array_deleter<OptimizerFunctionPtr>()), 0, 2, leftBound, rightBound);
	ags.SetTask(task);

	steady_clock::time_point start = steady_clock::now();
	for (int i = 1; i <= 100; i++)
	{
		function->SetFunctionNumber(i);
		globalMinPoint[0] = function->GetMinXCoordinate();
		globalMinPoint[1] = function->GetMinYCoordinate();

		auto result = ags.StartOptimization(globalMinPoint, optimizercore::StopCriterionType::Precision);
		auto stat = result.GetSolution();
		y = stat.GetOptimumPoint().get();
		double helderConst = ags.GetLipschitzConst(result.GetNumberOfFunctionals() - 1);
		avgHelderConst += helderConst / 100.0;

		printf("Function #%i\n", i);
		printf("Evaluated optimum point (%f  ; %f)\n", y[0], y[1]);
		printf("It_count: %i\n", result.GetNumberOfCalculations(result.GetNumberOfFunctionals() - 1));
		printf("Function value %f\n", stat.GetOptimumValue());
		printf("Helder const evaluation: %f", helderConst);
		printf("\n--------------------------------\n");
		meanItCount += result.GetNumberOfCalculations(result.GetNumberOfFunctionals() - 1) / 100.0;
		if (optimizercore::utils::NormNDimMax(y, globalMinPoint, 2) > 0.01)
		{
			err_count++;
			results[i - 1] = 21000000;
			//break;
		}
		else {
			results[i - 1] = result.GetNumberOfCalculations(result.GetNumberOfFunctionals() - 1);
			if (max_count < results[i - 1])
				max_count = results[i - 1];
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

		sprintf_s(filename, "VAGr R= %.3f map_type= %i local percent= %i, threadsNum= %i maps_num=%i.csv",
			*algParams.r, m_type, algParams.localMixParameter, algParams.numberOfThreads,
			algParams.numberOfMaps);

		printf("\nVAGr R= %f map_type= %i local percent= %i, threadsNum= %i\n",
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
