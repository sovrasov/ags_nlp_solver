#include <cmath>
#include <cassert>
#include <stdexcept>
#include "OptimizerAlgorithm.hpp"
#include "Map.h"
#include "CoreUtils.hpp"

using namespace optimizercore;
using namespace optimizercore::utils;

OptimizerAlgorithm::OptimizerAlgorithm() : r(nullptr)
{
	mCurrentStorageSize = 2000;
	mStorageReallocStep = 2000;
	mStorageReallocBarrier = 200;
	mIsAlgorithmMemoryAllocated = false;

	mLocalStartIterationNumber = 1;
	mNumberOfThreads = 1;
	mMaxNumberOfIterations = 5000;
	mNextPoints = nullptr;
	mRestrictions = nullptr;
	next_points = nullptr;
	mIntervalsForTrials = nullptr;

	r_point = 1; l_point = 0;

	mIsTaskInitialized = false;
	mIsParamsInitialized = false;
}

void OptimizerAlgorithm::SetTask(OptimizerTask task)
{
	mTask = task;
	
	mRestrictionsNumber = mTask.GetNumberOfRestrictions();
	mTargetFunction = mTask.GetTaskFunctions().get()[mRestrictionsNumber].get();

	if (mRestrictionsNumber > 0)
	{
		if (mRestrictions)
			delete[] mRestrictions;
		mRestrictions = new OptimizerFunction*[mRestrictionsNumber];
		for (int i = 0; i < mRestrictionsNumber; i++)
			mRestrictions[i] = mTask.GetTaskFunctions().get()[i].get();
	}

	mIsTaskInitialized = true;
}
OptimizerSearchSequence OptimizerAlgorithm::GetSearchSequence() const
{
	return OptimizerSearchSequence(sq, it_count + 1, dim,
		static_cast<MapType> (m_type), mMapTightness);
}
double OptimizerAlgorithm::GetLipschitzConst(int fNumber) const
{
	assert(fNumber >= 0 && fNumber <= restr_count);
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
		local_percent = params.localMixParameter;
		mLocalMixType = true;
	}
	else	{
		local_percent = 20 - params.localMixParameter;
		mLocalMixType = false;
	}
	alpha = params.localExponent;
	dim = params.algDimention;
	mMapTightness = params.mapTightness;
	m_type = static_cast<int>(params.mapType);
	mMaxNumberOfIterations = params.maxIterationsNumber;
	r = params.r;
	if (mNextPoints)
		utils::DeleteMatrix(mNextPoints, mNumberOfThreads);
	mNextPoints = utils::AllocateMatrix<double>(mNumberOfThreads, dim);
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

	if (m_type == 3)
		mStorageReallocBarrier = mNumberOfThreads*(int)std::pow(2, dim) + 100;
	else
		mStorageReallocBarrier = mNumberOfThreads + 10;

	sq[0].x = l_point;
	sq[1].x = r_point;
	sq[0].v = sq[1].v = 0;
	sq[0].val = sq[1].val = 0;
}
OptimizerResult OptimizerAlgorithm::StartOptimization(const double* a, StopCriterionType stopType)
{
	assert(mIsParamsInitialized && mIsTaskInitialized);
	assert(dim == mTask.GetTaskDimention());

	InitializeInformationStorages();
	
	v_max = 1;
	double preimages[32], val;
	bool stop = false;
	int p_count = 0, i_count=0, v = 1, currentThrNum = 1;
	
	next_points[0].x = (r_point + l_point) / 2;
	y = mNextPoints[0];
	mapd(next_points[0].x, mMapTightness, y, dim, m_type);

	for (it_count = 2; i_count < mMaxNumberOfIterations && !stop; it_count++)
	{
#pragma omp parallel for num_threads(currentThrNum)
		for (int i = 0; i < currentThrNum; i++)
		{
			next_points[i].v = GetIndex(next_points + i, mNextPoints[i]);
#pragma omp critical
			if (next_points[i].v > v)
				v = next_points[i].v;
		}
		if (m_type == 3)
		{
			i_count++;
			for (int i = 0; i < currentThrNum; i++)
			{
				y = mNextPoints[i];
				invmad(mMapTightness, preimages, 32, &p_count, y, dim, 1);
				val = next_points[i].val;
				for (int k = 0; k < p_count; k++)
				{
					next_points[i].x = preimages[k];
					ins_pos = insert(sq, next_points[i], it_count + k);
					UpdateLipConsts(v_indexes[v - 1], next_points[i].v);
					v_indexes[v - 1]->Add(sq[ins_pos]);
				}
				if (i == currentThrNum - 1)
					it_count--;
				it_count += p_count;
			}
		}
		else
		{
			i_count++;
			for (int i = 0; i < currentThrNum; i++)
			{
				ins_pos = insert(sq, next_points[i], it_count + i);
				UpdateLipConsts(v_indexes[sq[ins_pos].v - 1], sq[ins_pos].v);
				v_indexes[sq[ins_pos].v - 1]->Add(sq[ins_pos]);
			}
			it_count += currentThrNum - 1;
		}
		//////
//		if (v_indexes[v - 1]->GetSize() > 20)
//			r[v - 1] = 2;
		//////////

		if (v > v_max)
		{
			for (int i = 0; i < v - 1; i++)
				set_ranks[i] = -reserves*lip_const[i];
			v_max = v;
		}
		if (v == v_max)
			set_ranks[v - 1] = v_indexes[v - 1]->GetMinimumValue();

		if (i_count >= mLocalStartIterationNumber)	{
			if (i_count % (12 - local_percent) == 0 && local_percent > 0)
				UpdateRanks(mLocalMixType);
			else
				UpdateRanks(!mLocalMixType);
		}
		else
			UpdateRanks(false);

		if (i_count < mNumberOfThreads + 10)
			currentThrNum = 1;
		else
			currentThrNum = mNumberOfThreads;

		for (int i = 0; i < currentThrNum && !stop; i++)
		{
			OptimizerTrialPoint left = mIntervalsForTrials[i].left;
			OptimizerTrialPoint right = mIntervalsForTrials[i].right;

			if (left.v == right.v)
				next_points[i].x = (left.x + right.x) / 2
				- sgn(right.val - left.val)*pow((right.val - left.val)
				/ lip_const[right.v - 1], dim) / (2 * r[right.v - 1]);
			else
				next_points[i].x = (right.x + left.x) / 2;

			mapd(next_points[i].x, mMapTightness, mNextPoints[i], dim, m_type);
			y = mNextPoints[i];

			if (stopType == StopCriterionType::OptimalPoint)	{
				if (NormNDimMax(y, a, dim) < eps)
				{
					stop = true;
					result = next_points[i];
				}
			}
			else	{
				if (std::pow(right.x - left.x, 1.0 / dim) < eps)
				{
					stop = true;
					result = next_points[i];
				}
			}
		}

		if (it_count > mCurrentStorageSize - mStorageReallocBarrier)
		{
			reAllocMem(&sq, mCurrentStorageSize, mStorageReallocStep);
			mCurrentStorageSize += mStorageReallocStep;
		}
	}
	
	sq[it_count].x = result.x;
	sq[it_count].val = result.val = mTargetFunction->Calculate(y);
	sq[it_count].v = result.v;
	
	if (stopType == StopCriterionType::Precision)
	{
		if (result.v == mRestrictionsNumber + 1)
			v_indexes[0]->Add(sq[it_count]);
		result = v_indexes[mRestrictionsNumber]->GetMinimumPoint();
	}
	
	mapd(result.x, mMapTightness, y, dim, m_type);
	SharedVector optPoint(new double[dim], array_deleter<double>());
	memcpy(optPoint.get(), y, dim*sizeof(double));

	OptimizerSolution Solution(i_count, result.val, result.x, dim, optPoint);

	return OptimizerResult(Solution);
}
void OptimizerAlgorithm::SetThreadsNum(int num)
{
	if (num > 0 && num < 100)
	{
		if (mNextPoints != nullptr)
			utils::DeleteMatrix(mNextPoints, mNumberOfThreads);
		mNumberOfThreads = num;
		if (next_points)
			delete[] next_points;
		if (mIntervalsForTrials)
			delete[] mIntervalsForTrials;
		mIntervalsForTrials = new OptimaizerInterval[num];
		next_points = new OptimizerTrialPoint[num];
		mNextPoints = utils::AllocateMatrix<double>(mNumberOfThreads, dim);
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
	if (next_points)
		delete[] next_points;
	if (mIsAlgorithmMemoryAllocated)
	{
		delete[] sq;
		delete[] lip_const;
		delete[] set_ranks;
		for (int j = 0; j < mRestrictionsNumber + 1; j++)
			delete v_indexes[j];
		delete[] v_indexes;
	}
}
void OptimizerAlgorithm::UpdateLipConsts(IndxSet* set, int setNumber)
{
	double max = lip_const[setNumber - 1], tmp;
	int set_size = set->GetSize();
	OptimizerTrialPoint cur_point;

	for(int k = 0; k < set_size; k++)
	{
		cur_point = set->Get(k);
		if (sq[ins_pos].x != cur_point.x)
		{
			tmp = std::abs(sq[ins_pos].val - cur_point.val) 
				/ pow(std::abs(sq[ins_pos].x - cur_point.x), 1.0 / dim);
			if (max < tmp)
				max = tmp;
		}
	}
	if (max != 0)
		lip_const[setNumber - 1] = max;
	else
		lip_const[setNumber - 1] = 1;
}
void OptimizerAlgorithm::UpdateRanks(bool isLocal)
{
	double dx, curr_rank;

	for (int i = 0; i < mNumberOfThreads; i++)
		mIntervalsForTrials[i].rank = -HUGE_VAL;
	
	for (int j = 0; j < it_count; j++)
	{
		dx = pow(std::abs(sq[j + 1].x - sq[j].x), 1.0 / dim);

		if (dx == 0)
			throw std::runtime_error("Interval already exists!");

		if (sq[j + 1].v == sq[j].v)
		{
			int v = sq[j].v - 1;
			curr_rank = dx + Pow2((sq[j + 1].val - sq[j].val) / (r[v] * lip_const[v])) / dx
				- 2 * (sq[j + 1].val + sq[j].val - 2 * set_ranks[v]) / (r[v] * lip_const[v]);
		}
		else if (sq[j + 1].v > sq[j].v)
		{
			int v = sq[j + 1].v - 1;
			curr_rank = 2 * dx - 4 * (sq[j + 1].val - set_ranks[v]) / (r[v] * lip_const[v]);
		}
		else //if (sq[j + 1].v < sq[j].v)
		{
			int v = sq[j].v - 1;
			curr_rank = 2 * dx - 4 * (sq[j].val - set_ranks[sq[j].v - 1]) / (r[v] * lip_const[v]);
		}

		if (isLocal)
		{
			if (sq[j + 1].v == sq[j].v)
			{
				int v = sq[j].v - 1;
				curr_rank /= sqrt((sq[j + 1].val - set_ranks[v])*
					(sq[j].val - set_ranks[v])) / lip_const[v] + pow(1.5, -alpha);
			}
			else if (sq[j + 1].v > sq[j].v)
			{
				int v = sq[j + 1].v - 1;
				curr_rank /= (sq[j + 1].val - set_ranks[v])
					/ lip_const[v] + pow(1.5, -alpha);
			}
			else //if (sq[j + 1].v < sq[j].v)
			{
				int v = sq[j].v - 1;
				curr_rank /= (sq[j].val - set_ranks[v])
					/ lip_const[v] + pow(1.5, -alpha);
			}
		}
		
		if (curr_rank > mIntervalsForTrials[mNumberOfThreads - 1].rank)
		{
			OptimaizerInterval newInterval(sq[j], sq[j + 1], curr_rank);
			for (int i = 0; i < mNumberOfThreads; i++)
				if (mIntervalsForTrials[i].rank < newInterval.rank)
					std::swap(mIntervalsForTrials[i], newInterval);
		}
	}
}
int OptimizerAlgorithm::GetIndex(OptimizerTrialPoint* oneDimPoint, double* point)
{
	int indx = 1;
	for (int j = 0; j < mRestrictionsNumber; j++)
	{
		oneDimPoint->val = mRestrictions[j]->Calculate(point);
		if (oneDimPoint->val > 0)
			break;
		else
			indx++;
	}
	if (indx == mRestrictionsNumber + 1)
		oneDimPoint->val = mTargetFunction->Calculate(point);
	return indx;
}
void OptimizerAlgorithm::AllocMem()
{
	sq = new OptimizerTrialPoint[mCurrentStorageSize];
	//ranks = new double[mCurrentStorageSize];
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