#include <cmath>
#include <cassert>
#include <algorithm>

#include "OptimizerAlgorithm.hpp"
#include "CoreUtils.hpp"
#include "HookeJeevesLocalMethod.hpp"

using namespace optimizercore;
using namespace optimizercore::utils;

OptimizerAlgorithm::OptimizerAlgorithm()
{
	mCurrentStorageSize = 5000;
	mStorageReallocStep = 3000;
	mIsAlgorithmMemoryAllocated = false;

	mLocalStartIterationNumber = 1;
	mNumberOfThreads = 1;
	mMaxNumberOfIterations = 10000;
	mNextPoints = nullptr;
	mRestrictions = nullptr;
	mNextTrialsPoints = nullptr;
	mIntervalsForTrials = nullptr;
	r = nullptr;
	mPMap = nullptr;

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
	return OptimizerSearchSequence(mSearchInformationStorage, mMethodDimension,
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
	mMethodDimension = params.algDimention;
	mMapTightness = params.mapTightness;
	mMapType = static_cast<int>(params.mapType);

	if (mPMap)
		delete mPMap;
	if (mMapType < 4) {
		mPMap = new OptimizerMap(mMethodDimension, mMapTightness, params.mapType);
		mNumberOfMaps = 1;
	}
	else if (mMapType == 4)
		mPMap = new OptimizerMultiMap(MultimapType::Set, mMethodDimension, mMapTightness, params.numberOfMaps);
	else if(mMapType == 5)
		mPMap = new OptimizerMultiMap(MultimapType::Rotated, mMethodDimension, mMapTightness, params.numberOfMaps);
	mNumberOfMaps = params.numberOfMaps;
	mMaxNumberOfIterations = params.maxIterationsNumber;
	r = params.r;
	if (mNextPoints)
		utils::DeleteMatrix(mNextPoints, mNumberOfThreads);
	mNextPoints = utils::AllocateMatrix<double>(mNumberOfThreads, mMethodDimension);
	this->SetThreadsNum(params.numberOfThreads);

	mIsParamsInitialized = true;
}

void OptimizerAlgorithm::InitializeInformationStorages()
{
	if (mMapType == 4 || !mSpaceTransform.IsZeroConstraintActive()) {
		//insert the zero restriction for set map
		mRestrictionsNumber = mTask.GetNumberOfRestrictions() + 1;
		delete[] mRestrictions;
		mRestrictions = new OptimizerFunction*[mRestrictionsNumber];
		mRestrictions[0] = mSpaceTransform.GetZeroConstraint().get();
		for (int i = 1; i < mRestrictionsNumber; i++)
			mRestrictions[i] = mTask.GetTaskFunctions().get()[i - 1].get();
	}

	mFunctionalsCalculationsStat = new int[mRestrictionsNumber + 1];
	std::fill_n(mFunctionalsCalculationsStat, mRestrictionsNumber + 1, 0);

	if (!mIsAlgorithmMemoryAllocated) {
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

	if (mMapType > 3) {
		for (int i = 2; i <= mNumberOfMaps; i++)
			mSearchInformationStorage.emplace(i, 0.0, -1);
	}

	mNeedQueueRefill = true;
}

bool OptimizerAlgorithm::InsertNewTrials(int trailsNumber)
{
	bool storageInsertionError;
	
	if (mNumberOfMaps != 1 || mMapType == 3) {
		for (int i = 0; i < trailsNumber; i++) {
			if (mMapType == 5 || mMapType == 3 || mNextTrialsPoints[i].v != 0) {
				double preimages[MAX_PREIMAGES];
				int preimNumber = mPMap->GetAllPreimages(mNextPoints[i], preimages);
				for (int k = 0; k < preimNumber; k++) {
					mNextTrialsPoints[i].x = preimages[k];
					storageInsertionError =
						mSearchInformationStorage.insert(mNextTrialsPoints[i]).second;
					UpdateLipConsts(v_indexes[mNextTrialsPoints[i].v], mNextTrialsPoints[i]);
				}
			}
			else {
				storageInsertionError =
					mSearchInformationStorage.insert(mNextTrialsPoints[i]).second;
				UpdateLipConsts(v_indexes[mNextTrialsPoints[i].v], mNextTrialsPoints[i]);
			}
		}
	}
	else {
		for (int i = 0; i < trailsNumber; i++) {
			storageInsertionError =
				mSearchInformationStorage.insert(mNextTrialsPoints[i]).second;
			UpdateLipConsts(v_indexes[mNextTrialsPoints[i].v], mNextTrialsPoints[i]);
		}
	}

	return storageInsertionError;
}

OptimizerResult OptimizerAlgorithm::StartOptimization(
	const double* a, StopCriterionType stopType)
{
	assert(mIsParamsInitialized && mIsTaskInitialized);
	assert(mMethodDimension == mTask.GetTaskDimention());

	InitializeInformationStorages();
	
	double *y;
	bool stop = false;
	int iterationsCount = 0, v = 0, v_max = 0, 
		currentThrNum = 1, ranksUpdateErrCode;
	
	mNextTrialsPoints[0].x = 0.5;
	if (mMapType == 4)
		mNextTrialsPoints[0].x = pow(2, -(mMethodDimension + 1));

	mPMap->GetImage(mNextTrialsPoints[0].x, mNextPoints[0]);
	mSpaceTransform.Transform(mNextPoints[0], mNextPoints[0]);

	while (iterationsCount < mMaxNumberOfIterations && !stop)	{

		if (mMapType == 3 || mNumberOfMaps != 1 || mLocalMixParameter != 0)
			mNeedQueueRefill = true;
		iterationsCount++;

#pragma omp parallel for num_threads(currentThrNum)
		for (int i = 0; i < currentThrNum; i++)
		{
			mNextTrialsPoints[i].v = 
				GetIndex(mNextTrialsPoints + i, mNextPoints[i]);

			if (mMapType >= 3)
				mSpaceTransform.InvertTransform(mNextPoints[i], mNextPoints[i]);

#pragma omp critical
			{
				mFunctionalsCalculationsStat[mNextTrialsPoints[i].v]++;
				if (mNextTrialsPoints[i].v > v)
					v = mNextTrialsPoints[i].v;
			}
		}

		if (!InsertNewTrials(currentThrNum))
			break;

		if (v > v_max)
		{
			for (int i = 0; i < v; i++)
				set_ranks[i] = -reserves*lip_const[i];
			v_max = v;
			mNeedQueueRefill = true;
		}
		if (v == v_max)
			set_ranks[v] = v_indexes[v]->GetMinimumValue();

		//full update ranks or insert new intervals in queue
		if (mNeedQueueRefill) {
			if (iterationsCount >= mLocalStartIterationNumber) {
				if (iterationsCount % (12 - mLocalMixParameter) == 0
					&& mLocalMixParameter > 0)
					ranksUpdateErrCode = UpdateRanks(mLocalMixType);
				else
					ranksUpdateErrCode = UpdateRanks(!mLocalMixType);
			}
			else
				ranksUpdateErrCode = UpdateRanks(false);
		}
		else {
			for (int i = 0; i < currentThrNum; i++) {
				OptimizerInterval i1(mIntervalsForTrials[i].left, mNextTrialsPoints[i], 0, 0);
				OptimizerInterval i2(mNextTrialsPoints[i], mIntervalsForTrials[i].right, 0, 0);
				i1.R = CalculateGlobalR(i1);
				i2.R = CalculateGlobalR(i2);
				mQueue.Push(i1);
				mQueue.Push(i2);
			}
		}

		if (iterationsCount >= mNumberOfThreads + 10)
			currentThrNum = mNumberOfThreads;

		for (int i = 0; i < currentThrNum && !stop; i++)	{
			//get new interval from queue
			mIntervalsForTrials[i] = mQueue.Pop();
			OptimizerTrialPoint left = mIntervalsForTrials[i].left;
			OptimizerTrialPoint right = mIntervalsForTrials[i].right;

			//fill mIntervalsForTrials
			if (left.v == right.v)
				mNextTrialsPoints[i].x = (left.x + right.x) / 2
				- sgn(right.val - left.val)*pow(fabs(right.val - left.val)
				/ lip_const[right.v], mMethodDimension) / (2 * r[right.v]);
			else
				mNextTrialsPoints[i].x = (right.x + left.x) / 2;
			
			mPMap->GetImage(mNextTrialsPoints[i].x, mNextPoints[i]);
			mSpaceTransform.Transform(mNextPoints[i], mNextPoints[i]);
			
			y = mNextPoints[i];

			if (stopType == StopCriterionType::OptimalPoint)	{
				if (NormNDimMax(y, a, mMethodDimension) < eps)	{
					stop = true;
					mOptimumEvaluation = mNextTrialsPoints[i];
				}
			}
			else	{
				if (pow(right.x - left.x, 1.0 / mMethodDimension) < eps)	{
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
	
	mPMap->GetImage(mOptimumEvaluation.x, y);
	mSpaceTransform.Transform(y, y);
	
	SharedVector optPoint(new double[mMethodDimension], array_deleter<double>());
	std::copy_n(y, mMethodDimension, optPoint.get());

	OptimizerSolution solution(iterationsCount, mOptimumEvaluation.val,
		mOptimumEvaluation.x, mMethodDimension, optPoint);

	if (mNeedLocalVerification)
		return OptimizerResult(DoLocalVerification(solution),
			SharedIntVector(mFunctionalsCalculationsStat, array_deleter<int>()), mRestrictionsNumber + 1);
	else
		return OptimizerResult(solution,
			SharedIntVector(mFunctionalsCalculationsStat, array_deleter<int>()), mRestrictionsNumber + 1);
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
		bestLocalValue, 0.5, mMethodDimension, localOptimum);

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
		mIntervalsForTrials = new OptimizerInterval[num];
		mNextTrialsPoints = new OptimizerTrialPoint[num];
		mNextPoints = utils::AllocateMatrix<double>(
			mNumberOfThreads, mMethodDimension);
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
	if (mPMap)
		delete mPMap;
	if (mIsAlgorithmMemoryAllocated)
	{
		delete[] lip_const;
		delete[] set_ranks;
		for (int j = 0; j < mRestrictionsNumber + 1; j++)
			delete v_indexes[j];
		delete[] v_indexes;
	}
}
void OptimizerAlgorithm::UpdateLipConsts(MultimapIndxSet* set,
	const OptimizerTrialPoint& value)
{
	double max = lip_const[value.v];
	if (max == 1) max = 0;
	int setNumber = (int)value.x;
	int set_size = set->GetSize(setNumber);
	OptimizerTrialPoint cur_point;

	for (int k = 0; k < set_size; k++)
	{
		cur_point = set->Get(k, setNumber);
		if (value.x != cur_point.x)
			max = fmax(max,
			fabs(value.val - cur_point.val)
			/ pow(fabs(value.x - cur_point.x), 1.0 / mMethodDimension));
	}

	if (max != 0) {
		if (lip_const[value.v] != max)
			mNeedQueueRefill = true;
		lip_const[value.v] = max;
	}
	else
		lip_const[value.v] = 1;

	if (set->GetMinimumValue() > value.val)
		mNeedQueueRefill = true;

	set->Add(value);
}
double optimizercore::OptimizerAlgorithm::CalculateGlobalR(const OptimizerInterval& interval)
{
	double R;
	const OptimizerTrialPoint& left = interval.left;
	const OptimizerTrialPoint& right = interval.right;

	double dx = pow(right.x - left.x, 1.0 / mMethodDimension);

	if (dx == 0)
		return HUGE_VAL;

	if (right.v == left.v)	{
		int v = left.v;
		R = dx + Pow2((right.val - left.val)
			/ (r[v] * lip_const[v])) / dx
			- 2 * (right.val + left.val - 2 * set_ranks[v])
			/ (r[v] * lip_const[v]);
	}
	else if (right.v > left.v)	{
		int v = right.v;
		R = 2 * dx - 4 * (right.val - set_ranks[v]) /
			(r[v] * lip_const[v]);
	}
	else {//if (rightIt->v < leftIt->v)	
		int v = left.v;
		R = 2 * dx - 4 * (left.val - set_ranks[v])
			/ (r[v] * lip_const[v]);
	}

	return R;
}
int OptimizerAlgorithm::UpdateRanks(bool isLocal)
{
	double dx, curr_rank;

	mQueue.Clear();

	for (int i = 0; i < mNumberOfThreads; i++)
		mIntervalsForTrials[i].R = -HUGE_VAL;

	auto leftIt = mSearchInformationStorage.begin();
	auto rightIt = mSearchInformationStorage.begin();
	++rightIt;

	int storageSize = mSearchInformationStorage.size();
	
	for (int j = 0; j < storageSize - 1; j++)
	{
		dx = pow(rightIt->x - leftIt->x, 1.0 / mMethodDimension);

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
		
		mQueue.Push(OptimizerInterval(OptimizerTrialPoint(*leftIt),
			OptimizerTrialPoint(*rightIt), curr_rank, 0));
		/*
		if (curr_rank > mIntervalsForTrials[mNumberOfThreads - 1].R)
		{
			OptimizerInterval newInterval(
				OptimizerTrialPoint(*leftIt),
				OptimizerTrialPoint(*rightIt), curr_rank, 0);
			for (int i = 0; i < mNumberOfThreads; i++)
				if (mIntervalsForTrials[i].R < newInterval.R)
					std::swap(mIntervalsForTrials[i], newInterval);
		}
		*/
		++leftIt;
		++rightIt;
	}
	mNeedQueueRefill = false;
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
	v_indexes = new MultimapIndxSet*[mRestrictionsNumber + 1];

	for (int j = 0; j < mRestrictionsNumber + 1; j++)
	{
		set_ranks[j] = 0;
		lip_const[j] = 1;
		v_indexes[j] = new MultimapIndxSet(mCurrentStorageSize, mNumberOfMaps, mStorageReallocStep);
	}
}