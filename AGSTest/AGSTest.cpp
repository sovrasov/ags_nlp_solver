#include "OptimizerAlgorithm.hpp"
#include "CoreUtils.hpp"
#include <cmath>
#include "SequenceVisualizer.h"
#include "GKLSClassTest.h"
#include "TestTaskFactory.hpp"
#include "OptimizerSTLFunctionWrapper.hpp"
#include "HookeJeevesLocalMethod.hpp"

using namespace optimizercore;

int main(int argc, char* argv[])
{
	int map_type = 1, local_percent = 3, alpha = 15;

	double r = 4.5, eps = 0.001, res = 0;
	/*
	if (argc == 13)
	for (int i = 1; i<argc; i++)	{
	if (!strcmp(argv[i], "-r"))
	r = atof(argv[i + 1]);
	else if (!strcmp(argv[i], "-mt"))
	map_type = atoi(argv[i + 1]);
	else if (!strcmp(argv[i], "-p"))
	eps = atof(argv[i + 1]);
	else if (!strcmp(argv[i], "-lc"))
	local_percent = atof(argv[i + 1]);
	else if (!strcmp(argv[i], "-dim"))
	GKLS_dim = atoi(argv[i + 1]);
	else if (!strcmp(argv[i], "-alph"))
	alpha = atof(argv[i + 1]);
	}
	else
	{
	printf("\nUsage:\n-r reliability\n-mt map type\n-p precision\n-lc percent of local iterations\n-dim dimention of tasks\n-alph local parameter\n");
	system("PAUSE");
	exit(0);
	}
	*/

	int localStartIteration = 10;
	local_percent = 0;
	r = 3;
	eps = 0.01;
	map_type = 1;
	alpha = 20;
	int numOfThreads = 1;
	int taskdim = 2;
	gklsfunction::GKLSClass gklsClass = gklsfunction::GKLSClass::Simple;
	int maxIterNumber = 2000;
	//	ParseArguments(argc, argv, numOfThreads, taskdim, gklsClass, maxIterNumber);

	OptimizerAlgorithm ags;
	OptimizerTask task = TestTaskFactory::GetTask(3);

	taskdim = task.GetTaskDimention();
	OptimizerParameters params(maxIterNumber, numOfThreads, eps, &r, &res, taskdim, alpha, local_percent,
		localStartIteration,
		static_cast<MapType>(map_type), 12, false, LocalTuningMode::Adaptive);
	params.reserves = &res;
	params.r = new double[task.GetNumberOfRestrictions() + 1];
	for (int i = 0; i < task.GetNumberOfRestrictions() + 1; i++)
		params.r[i] = r;

	ags.SetParameters(params);
	ags.SetTask(task);

	double x_opt[5] = { -0.068, 1.962, 2.431, 9.833, 9.833 };

	auto result = ags.StartOptimization(x_opt, optimizercore::StopCriterionType::Precision);
	auto optPoint = result.GetSolution().GetOptimumPoint();

	/*
	localoptimizer::HookeJeevesLocalMethod localMethod;
	localMethod.SetEps(eps / 1000);
	localMethod.SetInitialStep(2*eps);
	localMethod.SetProblem(task);
	localMethod.SetStepMultiplier(2);
	localMethod.SetStartPoint(x_opt, task.GetTaskDimention());
	
	localMethod.StartOptimization(optPoint.get());
	*/
	for (int i = 0; i < params.algDimention; i++)
		printf("x[%i]: %f   ", i, optPoint.get()[i]);
	printf("\nFvalue %f\n", task.GetTaskFunctions().get()
		[task.GetNumberOfRestrictions()].get()->Calculate(x_opt));
	printf("iterations %i\n", result.GetSolution().GetIterationsCount());
	for (int i = 0; i <= task.GetNumberOfRestrictions(); i++)
	{
		printf("Evaluated helder constant #%i: %f\n", i + 1, ags.GetLipschitzConst(i));
		printf("Calculations counter for function #%i: %i\n", i + 1,
			((OptimizerSTLFunctionWrapper*)task.GetTaskFunctions().get()[i].get())->GetCalculationsCounter());
	}
	//VisualizeSolution(task, ags.GetSearchSequence(), result.GetSolution(), "st.png");
	TestGKLSClass(params, gklsClass, taskdim);

//	delete[] params.r;
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
	getchar();
	return 0;
}