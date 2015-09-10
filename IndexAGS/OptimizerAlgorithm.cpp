#include <cmath>
#include <cassert>
#include "OptimizerAlgorithm.hpp"
#include "Map.h"
#include "CoreUtils.hpp"
#include "HookeJeevesLocalMethod.hpp"

using namespace optimizercore;
using namespace optimizercore::utils;

OptimizerAlgorithm::OptimizerAlgorithm()
{
	mCurrentStorageSize = 2000;
	mStorageReallocStep = 2000;
	mIsAlgorithmMemoryAllocated = false;

	mLocalStartIterationNumber = 1;
	mNumberOfThreads = 1;
	mMaxNumberOfIterations = 5000;
	mNextPoints = nullptr;
	mRestrictions = nullptr;
	mNextTrialsPoints = nullptr;
	mIntervalsForTrials = nullptr;
	r = nullptr;

	mIsTaskInitialized = false;
	mIsParamsInitialized = false;
}

void OptimizerAlgorithm::SetTask(OptimizerTask task)
{
	mTask = task;
	
	mRestrictionsNumber = mTask.GetNumberOfRestrictions();
	mTargetFunction = mTask.GetTaskFunctions().get()[mRestrictionsNumber].get();
	mSpaceTransform = mTask.GetSpaceTransformation();

	OptimizerFunctionPtr zeroConstraint = mSpaceTransform.GetZeroConstraint();
	if (mSpaceTransform.IsZeroConstraintActive())
		mRestrictionsNumber++;

	if (mRestrictionsNumber > 0)
	{
		if (mRestrictions)
			delete[] mRestrictions;

		mRestrictions = new OptimizerFunction*[mRestrictionsNumber];
		if (mSpaceTransform.IsZeroConstraintActive())
		{
			mRestrictions[0] = zeroConstraint.get();
			for (int i = 1; i < mRestrictionsNumber; i++)
				mRestrictions[i] = mTask.GetTaskFunctions().get()[i - 1].get();
		}
		else
			for (int i = 0; i < mRestrictionsNumber; i++)
				mRestrictions[i] = mTask.GetTaskFunctions().get()[i].get();
	}

	mIsTaskInitialized = true;
}

OptimizerSearchSequence OptimizerAlgorithm::GetSearchSequence() const
{
	return OptimizerSearchSequence(mSearchInformationStorage, mMethodDimention,
		static_cast<MapType> (mMapType), mMapTightness, mSpaceTransform);
}

double OptimizerAlgorithm::GetLipschitzConst(int fNumber) const
{
	assert(fNumber >= 0 && fNumber <= mRestrictionsNumber);
	return lip_const[fNumber];
}

void OptimizerAlgorithm::SetParameters(OptimizerParameters params)
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

	mLocalStartIterationNumber = params.localAlgStartIterationNumber;
	reserves = params.reserves[0];
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
	r = params.r;
	if (mNextPoints)
		utils::DeleteMatrix(mNextPoints, mNumberOfThreads);
	mNextPoints = utils::AllocateMatrix<double>(mNumberOfThreads, mMethodDimention);
	this->SetThreadsNum(params.numberOfThreads);

	mIsParamsInitialized = true;
}

void OptimizerAlgorithm::InitializeInformationStorages()
{
	if (!mIsAlgorithmMemoryAllocated){
		AllocMem();
		mIsAlgorithmMemoryAllocated = true;
	}
	for (int j = 0; j < mRestrictionsNumber + 1; j++)
	{
		v_indexes[j]->Reset();
		set_ranks[j] = 0;
		lip_const[j] = 1;
	}

	mSearchInformationStorage.clear();
	mSearchInformationStorage.emplace(0.0, 0.0, -1);
	mSearchInformationStorage.emplace(1.0, 0.0, -1);
}

bool OptimizerAlgorithm::InsertNewTrials(int trailsNumber)
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
				storageInsertionError =
					mSearchInformationStorage.insert(mNextTrialsPoints[i]).second;
				UpdateLipConsts(
					v_indexes[mNextTrialsPoints[i].v], mNextTrialsPoints[i]);
			}
		}
	}
	else
		for (int i = 0; i < trailsNumber; i++)
		{
			storageInsertionError =
				mSearchInformationStorage.insert(mNextTrialsPoints[i]).second;
			UpdateLipConsts(v_indexes[mNextTrialsPoints[i].v], mNextTrialsPoints[i]);
		}
	return storageInsertionError;
}

OptimizerResult OptimizerAlgorithm::StartOptimization(
	const double* a, StopCriterionType stopType)
{
	assert(mIsParamsInitialized && mIsTaskInitialized);
	assert(mMethodDimention == mTask.GetTaskDimention());

	InitializeInformationStorages();
	
	double *y;
	bool stop = false;
	int iterationsCount = 0, v = 0, v_max = 0, 
		currentThrNum = 1, ranksUpdateErrCode;
	
	mNextTrialsPoints[0].x = 0.5;
	mapd(mNextTrialsPoints[0].x, mMapTightness, mNextPoints[0],
		mMethodDimention, mMapType);
	mSpaceTransform.Transform(mNextPoints[0], mNextPoints[0]);
	

	while (iterationsCount < mMaxNumberOfIterations && !stop)	{
		iterationsCount++;

#pragma omp parallel for num_threads(currentThrNum)
		for (int i = 0; i < currentThrNum; i++)
		{
			mNextTrialsPoints[i].v = 
				GetIndex(mNextTrialsPoints + i, mNextPoints[i]);

			if (mMapType == 3)
				mSpaceTransform.InvertTransform(mNextPoints[i], mNextPoints[i]);

#pragma omp critical
			if (mNextTrialsPoints[i].v > v)
				v = mNextTrialsPoints[i].v;
		}

		if (!InsertNewTrials(currentThrNum))
			break;

		if (v > v_max)
		{
			for (int i = 0; i < v; i++)
				set_ranks[i] = -reserves*lip_const[i];
			v_max = v;
		}
		if (v == v_max)
			set_ranks[v] = v_indexes[v]->GetMinimumValue();

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

			if (left.v == right.v)
				mNextTrialsPoints[i].x = (left.x + right.x) / 2
				- sgn(right.val - left.val)*pow(fabs(right.val - left.val)
				/ lip_const[right.v], mMethodDimention) / (2 * r[right.v]);
			else
				mNextTrialsPoints[i].x = (right.x + left.x) / 2;
			
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
	
	mOptimumEvaluation.v = GetIndex(&mOptimumEvaluation, y);
	mSearchInformationStorage.insert(mOptimumEvaluation);

	if (stopType == StopCriterionType::Precision)
	{
		if (mOptimumEvaluation.v == mRestrictionsNumber)
			v_indexes[mRestrictionsNumber]->Add(mOptimumEvaluation);
		mOptimumEvaluation = v_indexes[mRestrictionsNumber]->GetMinimumPoint();
	}
	
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
OptimizerSolution OptimizerAlgorithm::DoLocalVerification(OptimizerSolution startSolution)
{
	localoptimizer::HookeJeevesLocalMethod localMethod;
	localMethod.SetEps(eps / 100);
	localMethod.SetInitialStep(2*eps);
	localMethod.SetProblem(mTask);
	localMethod.SetStepMultiplier(2);
	localMethod.SetStartPoint(startSolution.GetOptimumPoint().get(),
		mTask.GetTaskDimention());

	SharedVector localOptimum(new double[mTask.GetTaskDimention()], array_deleter<double>());
	localMethod.StartOptimization(localOptimum.get());
	double bestLocalValue = mTargetFunction->Calculate(localOptimum.get());

	if (startSolution.GetOptimumValue() > bestLocalValue)
		return OptimizerSolution(startSolution.GetIterationsCount(),
		bestLocalValue, 0.5, mMethodDimention, localOptimum);

	return startSolution;
}
void OptimizerAlgorithm::SetThreadsNum(int num)
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
OptimizerAlgorithm::~OptimizerAlgorithm()
{
	if (mIntervalsForTrials)
		delete[] mIntervalsForTrials;
	if (mRestrictions)
		delete[] mRestrictions;
	if (mNextPoints)
		utils::DeleteMatrix(mNextPoints, mNumberOfThreads);
	if (mNextTrialsPoints)
		delete[] mNextTrialsPoints;
	if (mIsAlgorithmMemoryAllocated)
	{
		delete[] lip_const;
		delete[] set_ranks;
		for (int j = 0; j < mRestrictionsNumber + 1; j++)
			delete v_indexes[j];
		delete[] v_indexes;
	}
}
void OptimizerAlgorithm::UpdateLipConsts(IndxSet* set,
	const OptimizerTrialPoint& value)
{
	double max = lip_const[value.v];
	if (max == 1) max = 0;
	int set_size = set->GetSize();
	OptimizerTrialPoint cur_point;

	for (int k = 0; k < set_size; k++)
	{
		cur_point = set->Get(k);
		if (value.x != cur_point.x)
			max = fmax(max,
			fabs(value.val - cur_point.val)
			/ pow(fabs(value.x - cur_point.x), 1.0 / mMethodDimention));
	}

	if (max != 0)
		lip_const[value.v] = max;
	else
		lip_const[value.v] = 1;
	set->Add(value);
}
int OptimizerAlgorithm::UpdateRanks(bool isLocal)
{
	double dx, curr_rank;

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

		if (rightIt->v == leftIt->v)
		{
			int v = leftIt->v;
			curr_rank = dx + Pow2((rightIt->val - leftIt->val) 
				/ (r[v] * lip_const[v])) / dx
				- 2 * (rightIt->val + leftIt->val - 2 * set_ranks[v]) 
				/ (r[v] * lip_const[v]);
		}
		else if (rightIt->v > leftIt->v)
		{
			int v = rightIt->v;
			curr_rank = 2 * dx - 4 * (rightIt->val - set_ranks[v]) /
				(r[v] * lip_const[v]);
		}
		else //if (rightIt->v < leftIt->v)
		{
			int v = leftIt->v;
			curr_rank = 2 * dx - 4 * (leftIt->val - set_ranks[v]) 
				/ (r[v] * lip_const[v]);
		}

		if (isLocal)
		{
			if (rightIt->v == leftIt->v)
			{
				int v = rightIt->v;
				curr_rank /= sqrt((rightIt->val - set_ranks[v])*
					(leftIt->val - set_ranks[v])) / lip_const[v] + pow(1.5, -mAlpha);
			}
			else if (rightIt->v > leftIt->v)
			{
				int v = rightIt->v;
				curr_rank /= (rightIt->val - set_ranks[v])
					/ lip_const[v] + pow(1.5, -mAlpha);
			}
			else //if (rightIt->v < leftIt->v)
			{
				int v = leftIt->v;
				curr_rank /= (leftIt->val - set_ranks[v])
					/ lip_const[v] + pow(1.5, -mAlpha);
			}
		}

		if (curr_rank > mIntervalsForTrials[mNumberOfThreads - 1].rank)
		{
			OptimaizerInterval newInterval(
				OptimizerTrialPoint(*leftIt),
				OptimizerTrialPoint(*rightIt), curr_rank, 0);
			for (int i = 0; i < mNumberOfThreads; i++)
				if (mIntervalsForTrials[i].rank < newInterval.rank)
					std::swap(mIntervalsForTrials[i], newInterval);
		}
		++leftIt;
		++rightIt;
	}
	return 0;
}
int OptimizerAlgorithm::GetIndex(OptimizerTrialPoint* oneDimPoint, double* point)
{
	int indx = 0;
	for (int j = 0; j < mRestrictionsNumber; j++)
	{
		oneDimPoint->val = mRestrictions[j]->Calculate(point);
		if (oneDimPoint->val > 0)
			break;
		else
			indx++;
	}
	if (indx == mRestrictionsNumber)
		oneDimPoint->val = mTargetFunction->Calculate(point);
	return indx;
}
void OptimizerAlgorithm::AllocMem()
{
	lip_const = new double[mRestrictionsNumber + 1];
	set_ranks = new double[mRestrictionsNumber + 1];
	v_indexes = new IndxSet*[mRestrictionsNumber + 1];

	for (int j = 0; j < mRestrictionsNumber + 1; j++)
	{
		set_ranks[j] = 0;
		lip_const[j] = 1;
		v_indexes[j] = new IndxSet(mCurrentStorageSize, mStorageReallocStep);
	}
}