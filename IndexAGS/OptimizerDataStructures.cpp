#include <cmath>
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