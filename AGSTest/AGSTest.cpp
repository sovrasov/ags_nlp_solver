#include "OptimizerAlgorithm.hpp"
#include "CoreUtils.hpp"
#include <cmath>
#include "SequenceVisualizer.hpp"
#include "GKLSClassTest.hpp"
#include "TestTaskFactory.hpp"
#include "OptimizerSTLFunctionWrapper.hpp"
#include "HookeJeevesLocalMethod.hpp"
#include "TestGrishaginClass.hpp"
#include "HolderConstEvaluate.hpp"

using namespace optimizercore;

enum class TestSet {GKLS, Gris};

int main(int argc, char* argv[])
{
	int map_type = 1;
	int local_percent = 0;
	int alpha = 15;
	double r = 3.5, eps = 0.01, res = 0;
	int localStartIteration = 10;
	int numOfThreads = 1;
	int taskdim = 2;
	int maxIterNumber = 20000;
	int numberOfMaps = 2;

	gklsfunction::GKLSClass gklsClass = gklsfunction::GKLSClass::Simple;
	TestSet set = TestSet::Gris;

	for (int i = 1; i < argc; i++)	{
		if (!strcmp(argv[i], "-r"))
			r = atof(argv[i + 1]);
		else if (!strcmp(argv[i], "-mt"))
			map_type = atoi(argv[i + 1]);
		else if (!strcmp(argv[i], "-eps"))
			eps = atof(argv[i + 1]);
		else if(!strcmp(argv[i], "-nm"))
			numberOfMaps = atoi(argv[i + 1]);
		else if (!strcmp(argv[i], "-nt"))
			numOfThreads = atoi(argv[i + 1]);
		else if (!strcmp(argv[i], "-gkls"))
			set = TestSet::GKLS;
		else if (!strcmp(argv[i], "-gris"))
			TestSet set = TestSet::Gris;
	}

	OptimizerAlgorithm ags;
	OptimizerTask task = TestTaskFactory::GetTask(0);

	taskdim = task.GetTaskDimension();
	OptimizerParameters params(maxIterNumber, numOfThreads, eps, &r, &res, taskdim, alpha, local_percent,
		localStartIteration,
		static_cast<MapType>(map_type), 12, numberOfMaps, false, LocalTuningMode::Adaptive);
	params.reserves = &res;
	params.r = new double[task.GetNumberOfRestrictions() + 2];
	std::fill_n(params.r, task.GetNumberOfRestrictions() + 2, r);

	/*
	ags.SetParameters(params);
	ags.SetTask(task);

	double x_opt[5] = { -0.068, 1.962, 2.431, 9.833, 9.833 };

	auto result = ags.StartOptimization(x_opt, optimizercore::StopCriterionType::Precision);
	auto optPoint = result.GetSolution().GetOptimumPoint();

	localoptimizer::HookeJeevesLocalMethod localMethod;
	localMethod.SetEps(eps / 1000);
	localMethod.SetInitialStep(2*eps);
	localMethod.SetProblem(task);
	localMethod.SetStepMultiplier(2);
	localMethod.SetStartPoint(x_opt, task.GetTaskDimension());
	
	localMethod.StartOptimization(optPoint.get());
	for (int i = 0; i < params.algDimention; i++)
		printf("x[%i]: %f   ", i, optPoint.get()[i]);
	printf("\nFvalue %f\n", result.GetSolution().GetOptimumValue());
	printf("iterations %i\n", result.GetSolution().GetIterationsCount());
	for (int i = 0; i <= task.GetNumberOfRestrictions(); i++)
	{
		printf("Evaluated helder constant #%i: %f\n", i + 1, ags.GetLipschitzConst(i));
		printf("Calculations counter for function #%i: %i\n", i + 1,
			((OptimizerSTLFunctionWrapper*)task.GetTaskFunctions().get()[i].get())->GetCalculationsCounter());
	}
	*/

	if(set == TestSet::Gris)
		TestVAGrisClass(params);
	else
		TestMultimapsGKLSClass(params, gklsClass, 4);

	/*
	params.mapType = MapType::Rotated;
	params.numberOfMaps = 2;
	params.r[0] = 3.0;
	for (int i = 0; i <= 5; i++)
	{
		TestGKLSClass(params, gklsClass, 3);
		params.r[0] += 0.1;
	}
	params.numberOfMaps = 2;
	params.r[0] = 3.0;
	for (int i = 0; i <= 5; i++)
	{
		TestGKLSClass(params, gklsClass, 3);
		params.r[0] += 0.1;
	}
	params.numberOfMaps = 3;
	params.r[0] = 3.0;
	for (int i = 0; i <= 5; i++)
	{
		TestGKLSClass(params, gklsClass, 3);
		params.r[0] += 0.1;
	}
	params.mapType = MapType::Set;
	params.numberOfMaps = 2;
	params.r[0] = 3.0;
	for (int i = 0; i <= 5; i++)
	{
		TestGKLSClass(params, gklsClass, 3);
		params.r[0] += 0.1;
	}
	params.mapType = MapType::Set;
	params.numberOfMaps = 3;
	params.r[0] = 3.0;
	for (int i = 0; i <= 5; i++)
	{
		TestGKLSClass(params, gklsClass, 3);
		params.r[0] += 0.1;
	}
	*/
	delete[] params.r;

	//RunEvaluationExpOnGKLSClass(0, 0.00001);
	//RunEvaluationExpOnGrishaginClass(0, 0.00001);
	//VisualizeSolution(task, ags.GetSearchSequence(), result.GetSolution(), "st.png");

//	delete[] params.reserves;

	//////////////GLOBAL
	/*

	params.mapType = MapType::Noninjective;
	params.localMixParameter = 0;
	params.r[0] = 3.1;
	for (int i = 0; i <= 4; i++)
	{
		TestGKLSClass(params, gklsClass, taskdim);
		params.r[0] += 0.1;
	}
	
	params.localMixParameter = 7;

	params.r[0] = 3.1;
	for (int i = 0; i <= 4; i++)
	{
		TestGKLSClass(params, gklsClass, taskdim);
		params.r[0] += 0.1;
	}
	/////////////////LOCAL

	params.mapType = MapType::Simple;
	params.r[0] = 3.4;
	for (int i = 0; i <= 12; i++)
	{
		TestGKLSClass(params, gklsClass, taskdim);
		params.r[0] += 0.1;
	}
	
	params.mapType = MapType::Noninjective;
	params.r[0] = 3.4;
	for (int i = 0; i <= 12; i++)
	{
		TestGKLSClass(params, gklsClass, taskdim);
		params.r[0] += 0.1;
	}
	*/
/*
optimizercore::OptimizerAlgorithm ags;
int fNumber = 2;
double *x_opt = new double[2];
char fileName[128];
auto *function = new VAGRisFunctionWrapper();
function->SetFunctionNumber(fNumber);
x_opt[0] = function->GetMinXCoordinate() - 0.5;
x_opt[1] = function->GetMinYCoordinate() - 0.5;

auto parabaloidF = new FunctionPtrWrapper();
	parabaloidF->SetFunction(parabaloid);
	ags.SetThreadsNum(1);
	ags.SetParameters((OptimizerFunction**)&parabaloidF, 0, params);
	auto result1 = ags.StartOptimization(x_opt, optimizercore::StopCriterionType::Precision);
	auto sequence1 = ags.GetSearchSequence();
	sprintf(fileName, "1threads_%R=%f_EPS=%f.png", r, eps);
	VisualizeSequences(parabaloidF, sequence1, sequence1, fileName);

	ags.SetThreadsNum(1);
	ags.SetParameters((OptimizerFunction**)&function, 0, params);
	auto result1 = ags.StartOptimization(x_opt, optimizercore::StopCriterionType::OptimalPoint);
	auto sequence1 = ags.GetSearchSequence();

	ags.SetThreadsNum(2);
	auto result2 = ags.StartOptimization(x_opt, optimizercore::StopCriterionType::OptimalPoint);
	auto sequence2 = ags.GetSearchSequence();
	sprintf(fileName, "2threads_%R=%f_EPS=%f.png", r, eps);
	VisualizeSequences(function, sequence1, sequence2, fileName);

	ags.SetThreadsNum(3);
	result2 = ags.StartOptimization(x_opt, optimizercore::StopCriterionType::OptimalPoint);
	sequence2 = ags.GetSearchSequence();
	sprintf(fileName, "3threads_%R=%f_EPS=%f.png", r, eps);
	VisualizeSequences(function, sequence1, sequence2, fileName);

	ags.SetThreadsNum(4);
	result2 = ags.StartOptimization(x_opt, optimizercore::StopCriterionType::OptimalPoint);
	sequence2 = ags.GetSearchSequence();
	sprintf(fileName, "4threads_%R=%f_EPS=%f.png", r, eps);
	VisualizeSequences(function, sequence1, sequence2, fileName);
	*/
	//getchar();
	return 0;
}