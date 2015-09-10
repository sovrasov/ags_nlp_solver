#include "OptimizerAlgorithmUnconstrained.hpp"
#include "HookeJeevesLocalMethod.hpp"

#include <cassert>
#include <algorithm>

using namespace optimizercore;
using namespace optimizercore::utils;

OptimizerAlgorithmUnconstrained::OptimizerAlgorithmUnconstrained()
{
	mIsAlgorithmMemoryAllocated = false;

	mLocalStartIterationNumber = 1;
	mNumberOfThreads = 1;
	mMaxNumberOfIterations = 5000;
	mNextPoints = nullptr;
	mNextTrialsPoints = nullptr;
	mIntervalsForTrials = nullptr;
	r = 2;
	mLocalTuningMode = LocalTuningMode::None;

	mIsTaskInitialized = false;
	mIsParamsInitialized = false;
}

void OptimizerAlgorithmUnconstrained::SetTask(OptimizerFunctionPtr function,
	OptimizerSpaceTransformation spaceTransform)
{
	assert(function);

	mTargetFunctionSmartPtr = function;
	mTargetFunction = function.get();
	mSpaceTransform = spaceTransform;

	mIsTaskInitialized = true;
}

OptimizerSearchSequence OptimizerAlgorithmUnconstrained::GetSearchSequence() const
{
	return OptimizerSearchSequence(mSearchInformationStorage, mMethodDimention,
		static_cast<MapType> (mMapType), mMapTightness, mSpaceTransform);
}

double OptimizerAlgorithmUnconstrained::GetLipschitzConst() const
{
	return mGlobalM;
}

void OptimizerAlgorithmUnconstrained::SetParameters(OptimizerParameters params)
{
	assert(params.algDimention);
	assert(params.eps > 0);
	assert(params.localAlgStartIterationNumber > 0);
	assert(params.mapTightness > 5 && params.mapTightness <= 20);
	assert(params.maxIterationsNumber > 0);
	assert(params.localMixParameter >= 0 && params.localMixParameter <= 20);
	assert(params.r != nullptr);
	assert(params.numberOfThreads > 0);
	assert(params.reserves != nullptr);
	assert(params.adaptiveLocalTuningParameter >= 0 &&
		params.adaptiveLocalTuningParameter <= 1);

	mLocalStartIterationNumber = params.localAlgStartIterationNumber;
	eps = params.eps;
	if (params.localMixParameter <= 10)	{
		mLocalMixParameter = params.localMixParameter;
		mLocalMixType = true;
	}
	else	{
		mLocalMixParameter = 20 - params.localMixParameter;
		mLocalMixType = false;
	}
	mNeedLocalVerification = params.localVerification;
	mAlpha = params.localExponent;
	mMethodDimention = params.algDimention;
	mMapTightness = params.mapTightness;
	mMapType = static_cast<int>(params.mapType);
	mMaxNumberOfIterations = params.maxIterationsNumber;
	mLocalTuningMode = params.localTuningMode;
	r = *params.r;
	if (mNextPoints)
		utils::DeleteMatrix(mNextPoints, mNumberOfThreads);
	mNextPoints = utils::AllocateMatrix<double>(mNumberOfThreads, mMethodDimention);
	this->SetThreadsNum(params.numberOfThreads);

	mIsParamsInitialized = true;
}

void OptimizerAlgorithmUnconstrained::InitializeInformationStorage()
{
	if (!mIsAlgorithmMemoryAllocated){
		AllocMem();
		mIsAlgorithmMemoryAllocated = true;
	}

	mZ = HUGE_VAL;
	mGlobalM = 1;
	mMaxIntervalNorm = 0;

	mSearchInformationStorage.clear();

	mapd(0.0, mMapTightness, mNextPoints[0], mMethodDimention, mMapType);
	mSpaceTransform.Transform(mNextPoints[0], mNextPoints[0]);
	mSearchInformationStorage.emplace(0.0, mTargetFunction->Calculate(mNextPoints[0]), 0);

	mapd(1.0, mMapTightness, mNextPoints[0], mMethodDimention, mMapType);
	mSpaceTransform.Transform(mNextPoints[0], mNextPoints[0]);
	mSearchInformationStorage.emplace(1.0, mTargetFunction->Calculate(mNextPoints[0]), 0);
}

bool OptimizerAlgorithmUnconstrained::InsertNewTrials(int trailsNumber)
{
	bool storageInsertionError;
	if (mMapType == 3)
	{
		int preimagesNumber = 0;
		double preimages[32];
		for (int i = 0; i < trailsNumber; i++)
		{
			invmad(mMapTightness, preimages, 32,
				&preimagesNumber, mNextPoints[i], mMethodDimention, 4);
			for (int k = 0; k < preimagesNumber; k++)
			{
				mNextTrialsPoints[i].x = preimages[k];
				auto insertionResult =
					mSearchInformationStorage.insert(mNextTrialsPoints[i]);

				if (!(storageInsertionError = insertionResult.second))
					break;

				UpdateGlobalM(insertionResult.first);
			}
		}
	}
	else
		for (int i = 0; i < trailsNumber; i++)
		{
			auto insertionResult =
				mSearchInformationStorage.insert(mNextTrialsPoints[i]);

			if (!(storageInsertionError = insertionResult.second))
				break;
			
			UpdateGlobalM(insertionResult.first);
		}
	return storageInsertionError;
}

OptimizerResult OptimizerAlgorithmUnconstrained::StartOptimization(
	const double* a, StopCriterionType stopType)
{
	assert(mIsParamsInitialized && mIsTaskInitialized);
	assert(mSpaceTransform.GetDomainDimention() == mMethodDimention);

	InitializeInformationStorage();

	double *y;
	bool stop = false;
	int iterationsCount = 0,
		currentThrNum = 1, ranksUpdateErrCode;

	mNextTrialsPoints[0].x = 0.5;
	mapd(mNextTrialsPoints[0].x, mMapTightness, mNextPoints[0],
		mMethodDimention, mMapType);
	mSpaceTransform.Transform(mNextPoints[0], mNextPoints[0]);

	while (iterationsCount < mMaxNumberOfIterations && !stop)	{
		iterationsCount++;

#pragma omp parallel for num_threads(currentThrNum)
		for (int i = 0; i < currentThrNum; i++)	{
			mNextTrialsPoints[i].val = mTargetFunction->Calculate(mNextPoints[i]);
			if (mMapType == 3)
				mSpaceTransform.InvertTransform(mNextPoints[i], mNextPoints[i]);
#pragma omp critical
			if (mNextTrialsPoints[i].val < mZ)
				mZ = mNextTrialsPoints[i].val;
		}

		if (!InsertNewTrials(currentThrNum))
			break;

		if (iterationsCount >= mLocalStartIterationNumber)	{
			if (iterationsCount % (12 - mLocalMixParameter) == 0
				&& mLocalMixParameter > 0)
				ranksUpdateErrCode = UpdateRanks(mLocalMixType);
			else
				ranksUpdateErrCode = UpdateRanks(!mLocalMixType);
		}
		else
			ranksUpdateErrCode = UpdateRanks(false);

		if (iterationsCount >= mNumberOfThreads + 10)
			currentThrNum = mNumberOfThreads;

		for (int i = 0; i < currentThrNum && !stop; i++)	{
			OptimizerTrialPoint left = mIntervalsForTrials[i].left;
			OptimizerTrialPoint right = mIntervalsForTrials[i].right;

			mNextTrialsPoints[i].x = (left.x + right.x) / 2
				- sgn(right.val - left.val)*pow(fabs(right.val - left.val)
				/ mIntervalsForTrials[i].localM, mMethodDimention) / (2 * r);

			mapd(mNextTrialsPoints[i].x, mMapTightness, mNextPoints[i],
				mMethodDimention, mMapType);
			mSpaceTransform.Transform(mNextPoints[i], mNextPoints[i]);

			y = mNextPoints[i];

			if (stopType == StopCriterionType::OptimalPoint)	{
				if (NormNDimMax(y, a, mMethodDimention) < eps)	{
					stop = true;
					mOptimumEvaluation = mNextTrialsPoints[i];
				}
			}
			else	{
				if (pow(right.x - left.x, 1.0 / mMethodDimention) < eps)	{
					stop = true;
					mOptimumEvaluation = mNextTrialsPoints[i];
				}
			}
		}
	}

	mOptimumEvaluation.val = mTargetFunction->Calculate(y);
	mSearchInformationStorage.insert(mOptimumEvaluation);

	if (stopType == StopCriterionType::Precision)
		mOptimumEvaluation = *std::min_element(mSearchInformationStorage.begin(),
			mSearchInformationStorage.cend(),
			[](OptimizerTrialPoint p1, OptimizerTrialPoint p2)
		{
			return p1.val < p2.val;
		});

	mapd(mOptimumEvaluation.x, mMapTightness, y, mMethodDimention, mMapType);
	mSpaceTransform.Transform(y, y);

	SharedVector optPoint(new double[mMethodDimention], array_deleter<double>());
	std::memcpy(optPoint.get(), y, mMethodDimention*sizeof(double));

	OptimizerSolution solution(iterationsCount, mOptimumEvaluation.val,
		mOptimumEvaluation.x, mMethodDimention, optPoint);

	if (mNeedLocalVerification)
		return OptimizerResult(DoLocalVerification(solution));
	else
		return OptimizerResult(solution);
}
OptimizerSolution OptimizerAlgorithmUnconstrained::DoLocalVerification(OptimizerSolution startSolution)
{
	OptimizerFunctionPtr *functions = new OptimizerFunctionPtr[1];
	functions[0] = mTargetFunctionSmartPtr;

	OptimizerTask localTask(std::shared_ptr<OptimizerFunctionPtr>(functions,
		utils::array_deleter<OptimizerFunctionPtr>()),
		0, mMethodDimention, mSpaceTransform.GetLeftDomainBound(),
		mSpaceTransform.GetRightDomainBound());

	localoptimizer::HookeJeevesLocalMethod localMethod;
	localMethod.SetEps(eps / 100);
	localMethod.SetInitialStep(2 * eps);
	localMethod.SetProblem(localTask);
	localMethod.SetStepMultiplier(2);
	localMethod.SetStartPoint(startSolution.GetOptimumPoint().get(),
		localTask.GetTaskDimention());

	SharedVector localOptimum(new double[mMethodDimention], array_deleter<double>());
	localMethod.StartOptimization(localOptimum.get());
	double bestLocalValue = mTargetFunction->Calculate(localOptimum.get());

	if (startSolution.GetOptimumValue() > bestLocalValue)
		return OptimizerSolution(startSolution.GetIterationsCount(),
		bestLocalValue, 0.5, mMethodDimention, localOptimum);

	return startSolution;
}
void OptimizerAlgorithmUnconstrained::SetThreadsNum(int num)
{
	if (num > 0 && num < 100)
	{
		if (mNextPoints != nullptr)
			utils::DeleteMatrix(mNextPoints, mNumberOfThreads);
		mNumberOfThreads = num;
		if (mNextTrialsPoints)
			delete[] mNextTrialsPoints;
		if (mIntervalsForTrials)
			delete[] mIntervalsForTrials;
		mIntervalsForTrials = new OptimaizerInterval[num];
		mNextTrialsPoints = new OptimizerTrialPoint[num];
		mNextPoints = utils::AllocateMatrix<double>(
			mNumberOfThreads, mMethodDimention);
	}
}
OptimizerAlgorithmUnconstrained::~OptimizerAlgorithmUnconstrained()
{
	if (mIntervalsForTrials)
		delete[] mIntervalsForTrials;
	if (mNextPoints)
		utils::DeleteMatrix(mNextPoints, mNumberOfThreads);
	if (mNextTrialsPoints)
		delete[] mNextTrialsPoints;
	if (mIsAlgorithmMemoryAllocated)
	{
	}
}
void OptimizerAlgorithmUnconstrained::UpdateGlobalM(
	std::set<OptimizerTrialPoint>::iterator& newPointIt)
{
	double max = mGlobalM;
	if (max == 1) max = 0;

	
	auto leftPointIt = newPointIt;
	auto rightPointIt = newPointIt;
	--leftPointIt;
	++rightPointIt;

	double leftIntervalNorm = pow(newPointIt->x - leftPointIt->x, 1.0 / mMethodDimention);
	double rightIntervalNorm = pow(rightPointIt->x - newPointIt->x, 1.0 / mMethodDimention);


	max = fmax(fmax(fabs(newPointIt->val - leftPointIt->val) / leftIntervalNorm,
		fabs(rightPointIt->val - newPointIt->val) /	rightIntervalNorm), max);
	
	
	mMaxIntervalNorm = 0;
	auto currentPointIt = mSearchInformationStorage.begin();
	auto nextPointIt = currentPointIt;
	++nextPointIt;

	while (nextPointIt != mSearchInformationStorage.cend())
	{
	//	if (currentPointIt->x != newPointIt->x)
	//		max = fmax(fabs(newPointIt->val - currentPointIt->val) /
	//			pow(fabs(newPointIt->x - currentPointIt->x), 1.0 / mMethodDimention), max);

		if (mLocalTuningMode != LocalTuningMode::None)
		//if (mLocalTuningMode == LocalTuningMode::Maximum)
			mMaxIntervalNorm = fmax(
				pow(nextPointIt->x - currentPointIt->x, 1.0 / mMethodDimention),
				mMaxIntervalNorm);

		++currentPointIt;
		++nextPointIt;
	}
	//if (currentPointIt->x != newPointIt->x)
	//	max = fmax(fabs(newPointIt->val - currentPointIt->val) /
	//	pow(fabs(newPointIt->x - currentPointIt->x), 1.0 / mMethodDimention), max);
	
	//printf("%e |", mMaxIntervalNorm);
	if (max != 0)
		mGlobalM = max;
	else
		mGlobalM = 1;
}
int OptimizerAlgorithmUnconstrained::UpdateRanks(bool isLocal)
{
	double dx, curr_rank, mu1 = -HUGE_VAL, localM = mGlobalM;
	double localMConsts[3];

	for (int i = 0; i < mNumberOfThreads; i++)
		mIntervalsForTrials[i].rank = -HUGE_VAL;

	auto leftIt = mSearchInformationStorage.begin();
	auto rightIt = mSearchInformationStorage.begin();
	++rightIt;

	int storageSize = mSearchInformationStorage.size();

	for (int j = 0; j < storageSize - 1; j++)
	{
		dx = pow(rightIt->x - leftIt->x, 1.0 / mMethodDimention);

		if (dx == 0)
			return 1;

		if (mLocalTuningMode != LocalTuningMode::None)	{
			std::set<OptimizerTrialPoint>::iterator rightRightIt = rightIt;

			if (j > 0 && j < storageSize - 2)	{
				++rightRightIt;

				std::swap(localMConsts[0], localMConsts[1]);
				std::swap(localMConsts[1], localMConsts[2]);

				localMConsts[2] = fabs(rightRightIt->val - rightIt->val) 
					/ pow(rightRightIt->x - rightIt->x, 1.0 / mMethodDimention);

				mu1 = fmax(fmax(localMConsts[0], localMConsts[1]), localMConsts[2]);
			}
			else if (j == 0)	{
				++rightRightIt;

				localMConsts[1] = fabs(rightIt->val - leftIt->val) / dx;
				localMConsts[2] = fabs(rightRightIt->val - rightIt->val) /
					pow(rightRightIt->x - rightIt->x, 1.0 / mMethodDimention);
				mu1 = fmax(localMConsts[1], localMConsts[2]);
			}
			else
				mu1 = fmax(localMConsts[1], localMConsts[2]);

			double mu2 = mGlobalM*dx / mMaxIntervalNorm;

			if (mLocalTuningMode == LocalTuningMode::Maximum)	{
				localM = fmax(fmax(mu1, mu2), 0.01);
			}
			else// LocalTuningMode::Adaptive
				localM = fmax(mu1*(1 - dx / mMaxIntervalNorm) + mu2, 0.01);
			//localM = fmax(mu1*mMConvolution + (1 - mMConvolution)*mu2, 0.01);
			//localM = fmax(mu1 / r + (1 - 1 / r)*mGlobalM, 0.01);
			//if (mu1 < 0 || mu2 < 0)
			//	throw - 1;
			//printf(" %e  %e ||", mu1, mu2);
		}

		curr_rank = dx + Pow2((rightIt->val - leftIt->val) / (r * localM)) / dx
			- 2 * (rightIt->val + leftIt->val - 2 * mZ) / (r * localM);
		if (isLocal)
			curr_rank /= sqrt((rightIt->val - mZ)*
			(leftIt->val - mZ)) / localM + pow(1.5, -mAlpha);

		if (curr_rank > mIntervalsForTrials[mNumberOfThreads - 1].rank)
		{
			OptimaizerInterval newInterval(
				OptimizerTrialPoint(*leftIt),
				OptimizerTrialPoint(*rightIt), curr_rank, localM);
			for (int i = 0; i < mNumberOfThreads; i++)
				if (mIntervalsForTrials[i].rank < newInterval.rank)
					std::swap(mIntervalsForTrials[i], newInterval);
		}
		++leftIt;
		++rightIt;
	}
	return 0;
}
void OptimizerAlgorithmUnconstrained::AllocMem()
{
}
