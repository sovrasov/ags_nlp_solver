#include <cmath>
#include <algorithm>
#include "OptimizerDataStructures.hpp"
#include "CoreUtils.hpp"

using namespace optimizercore;
using namespace optimizercore::utils;

IndxSet::IndxSet() : mMem(nullptr), mCurrentZMin(HUGE_VAL), mCurrentSetSize(0) {}

IndxSet::IndxSet(int max_size, int mem_st) : mCurrentSetSize(0), 
											mMemSize(max_size), mMemReallocationStep(mem_st)
{
	mMem = new OptimizerTrialPoint[mMemSize];
}

void IndxSet::Reset()
{
	mCurrentSetSize = 0;
	mCurrentZMin = HUGE_VAL;
}

double IndxSet::GetMinimumValue() const
{
	return mCurrentZMin;
}

OptimizerTrialPoint IndxSet::GetMinimumPoint() const
{
	OptimizerTrialPoint min;
	min = mMem[0];
	for (int i = 1; i < mCurrentSetSize; i++)
		if (min.val > mMem[i].val)
			min = mMem[i];

	return min;
}

void IndxSet::Add(const OptimizerTrialPoint& trial)
{
	if (trial.val < mCurrentZMin)
		mCurrentZMin = trial.val;

	mMem[mCurrentSetSize++] = trial;

	if (mCurrentSetSize == mMemSize)
	{
		reAllocMem(&mMem, mMemSize, mMemReallocationStep);
		mMemSize += mMemReallocationStep;
	}
}

int IndxSet::GetSize() const
{
	return mCurrentSetSize;
}

IndxSet::~IndxSet()
{
	if (mMem != nullptr)
		delete[] mMem;
}

OptimizerTrialPoint& IndxSet::operator[](int i)
{
	return mMem[i];
}

OptimizerTrialPoint IndxSet::Get(int i) const
{
	return mMem[i];
}
//---------------------------------------------------------------------

optimizercore::MultimapIndxSet::MultimapIndxSet() : mMem(nullptr)
{
}

optimizercore::MultimapIndxSet::MultimapIndxSet(int memSize, int mapsNumber, int memReallocStep)
{
	mMem = new OptimizerTrialPoint*[mapsNumber];
	mMaxSetsSizes = new int[mapsNumber];
	mCurrentSetsSizes = new int[mapsNumber];

	for (int i = 0; i < mapsNumber; i++) {
		mMem[i] = new OptimizerTrialPoint[memSize];
		mMaxSetsSizes[i] = memSize;
		mCurrentSetsSizes[i] = 0;
	}

	mMemReallocationStep = memReallocStep;
	mNumberOfSubsets = mapsNumber;
	mCurrentZMin = HUGE_VAL;
}

void optimizercore::MultimapIndxSet::Add(const OptimizerTrialPoint & trial)
{
	int setNumber = (int)trial.x;
	mMem[setNumber][mCurrentSetsSizes[setNumber]++] = trial;

	if (trial.val < mCurrentZMin)
		mCurrentZMin = trial.val;

	if (mCurrentSetsSizes[setNumber] == mMaxSetsSizes[setNumber]) {
		reAllocMem(&mMem[setNumber], mCurrentSetsSizes[setNumber], mMemReallocationStep);
		mMaxSetsSizes[setNumber] += mMemReallocationStep;
	}
}

OptimizerTrialPoint optimizercore::MultimapIndxSet::Get(int i, int mapNumber) const
{
	return mMem[mapNumber][i];
}

int optimizercore::MultimapIndxSet::GetSize(int mapNumber) const
{
	return mCurrentSetsSizes[mapNumber];
}

double optimizercore::MultimapIndxSet::GetMinimumValue() const
{
	return mCurrentZMin;
}

OptimizerTrialPoint optimizercore::MultimapIndxSet::GetMinimumPoint() const
{
	OptimizerTrialPoint min(0, HUGE_VAL, 0);

	for (int i = 0; i < mNumberOfSubsets; i++)
		for (int j = 0; j < mCurrentSetsSizes[i]; j++)
			if (min.val > mMem[i][j].val)
				min = mMem[i][j];

	return min;
}

void optimizercore::MultimapIndxSet::Reset()
{
	std::fill_n(mCurrentSetsSizes, mNumberOfSubsets, 0);
	mCurrentZMin = HUGE_VAL;
}

optimizercore::MultimapIndxSet::~MultimapIndxSet()
{
	if (mMem != nullptr) {
		for (int i = 0; i < mNumberOfSubsets; i++)
			delete[] mMem[i];
		delete[] mMem;
		delete[] mCurrentSetsSizes;
		delete[] mMaxSetsSizes;
	}
}