#include "solver.hpp"
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>

namespace
{
    const double zeroHLevel = 1e-12;
}

NLPSolver::NLPSolver() {}

void NLPSolver::SetParameters(const SolverParameters& params)
{
  mParameters = params;
}

void NLPSolver::SetProblem(std::shared_ptr<IGOProblem<double>> problem)
{
  mProblem = problem;
}

std::vector<unsigned> NLPSolver::GetCalculationsStatistics() const
{
  return mCalculationsCounters;
}

std::vector<double> NLPSolver::GetHolderConstantsEstimations() const
{
  return mHEstimations;
}

void NLPSolver::InitDataStructures()
{
  double leftDomainBound[solverMaxDim], rightDomainBound[solverMaxDim];
  mProblem->GetBounds(leftDomainBound, rightDomainBound);
  mEvolvent = Evolvent(mProblem->GetDimension(), mParameters.evolventTightness,
    leftDomainBound, rightDomainBound);

  mNextPoints.resize(mParameters.numThreads);
  mOptimumEstimation.idx = -1;

  mZEstimations.resize(mProblem->GetConstraintsNumber() + 1);
  std::fill(mZEstimations.begin(), mZEstimations.end(),
            std::numeric_limits<double>::max());
  mNextIntervals.resize(mParameters.numThreads);
  mHEstimations.resize(mProblem->GetConstraintsNumber() + 1);
  std::fill(mHEstimations.begin(), mHEstimations.end(), 1.0);
  mCalculationsCounters.resize(mProblem->GetConstraintsNumber() + 1);
  std::fill(mCalculationsCounters.begin(), mCalculationsCounters.end(), 0);
  mQueue = PriorityQueue();
  mIterationsCounter = 0;
  mMinDelta = std::numeric_limits<double>::max();
  mMaxIdx = -1;
}

void NLPSolver::ClearDataStructures()
{
  for (const auto& ptr : mSearchInformation)
    delete ptr;
  mSearchInformation.clear();
  mQueue = PriorityQueue();
}

Trial NLPSolver::Solve()
{
  bool needStop = false;
  InitDataStructures();
  FirstIteration();

  do {
    EstimateOptimum();
    InsertIntervals();
    if (mNeedRefillQueue || mQueue.size() < mParameters.numThreads)
      RefillQueue();
    CalculateNextPoints();
    MakeTrials();
    needStop = mMinDelta < mParameters.eps;
    mIterationsCounter++;
  } while(mIterationsCounter < mParameters.trialsLimit && !needStop);

  ClearDataStructures();
  return mOptimumEstimation;
}

void NLPSolver::FirstIteration()
{
  Trial leftBound = Trial(0);
  leftBound.idx = -1;
  Trial rightBound = Trial(1.);
  rightBound.idx = -1;

  for (size_t i = 1; i <= mParameters.numThreads; i++)
  {
    mNextPoints[i - 1] = Trial((double)i / (mParameters.numThreads + 1));
    mEvolvent.GetImage(mNextPoints[i - 1].x, mNextPoints[i - 1].y);
  }

  MakeTrials();
  EstimateOptimum();

  for (size_t i = 0; i <= mParameters.numThreads; i++)
  {
    Interval* pNewInterval;
    if (i == 0)
      pNewInterval = new Interval(leftBound, mNextPoints[i]);
    else if (i == mParameters.numThreads)
      pNewInterval = new Interval(mNextPoints[i - 1], rightBound);
    else
      pNewInterval = new Interval(mNextPoints[i - 1], mNextPoints[i]);
    pNewInterval->delta = pow(pNewInterval->pr.x - pNewInterval->pl.x,
                              1. / mProblem->GetDimension());
    mMinDelta = std::min(mMinDelta, pNewInterval->delta);
    auto insRes = mSearchInformation.insert(pNewInterval);
    UpdateH(insRes.first);
  }
  RefillQueue();
  CalculateNextPoints();
  MakeTrials();
}

void NLPSolver::MakeTrials()
{
  for (size_t i = 0; i < mNextPoints.size(); i++)
  {
    int idx = 0;
    while(idx < mProblem->GetConstraintsNumber())
    {
      mNextPoints[i].idx = idx;
      double val = mProblem->Calculate(mNextPoints[i].y, idx);
      mCalculationsCounters[idx]++;
      if (val > 0)
      {
        mNextPoints[i].g[idx] = val;
        break;
      }
      idx++;
    }

    if(idx > mMaxIdx)
    {
      mMaxIdx = idx;
      for(int i = 0; i < mMaxIdx; i++)
        mZEstimations[i] = -mParameters.rEps;
      mNeedRefillQueue = true;
    }
    if (idx == mProblem->GetConstraintsNumber())
    {
      mCalculationsCounters[idx]++;
      mNextPoints[i].idx = idx;
      mNextPoints[i].g[idx] = mProblem->Calculate(mNextPoints[i].y, idx);
    }
    if(mNextPoints[i].idx == mMaxIdx &&
       mNextPoints[i].g[mMaxIdx] < mZEstimations[mMaxIdx])
    {
      mZEstimations[mMaxIdx] = mNextPoints[i].g[mMaxIdx];
      mNeedRefillQueue = true;
    }
  }
}

void NLPSolver::InsertIntervals()
{
  for (size_t i = 0; i < mParameters.numThreads; i++)
  {
    Interval* pOldInterval = mNextIntervals[i];
    Interval* pNewInterval = new Interval(mNextPoints[i], pOldInterval->pr);
    pOldInterval->pr = mNextPoints[i];
    pOldInterval->delta = pow(pOldInterval->pr.x - pOldInterval->pl.x,
                              1. / mProblem->GetDimension());
    pNewInterval->delta = pow(pNewInterval->pr.x - pNewInterval->pl.x,
                              1. / mProblem->GetDimension());
    mMinDelta = std::min(mMinDelta, pNewInterval->delta);
    mMinDelta = std::min(mMinDelta, pOldInterval->delta);

    auto insResult = mSearchInformation.insert(pNewInterval);
    bool wasInserted = insResult.second;
    if(!wasInserted)
      throw std::runtime_error("Error during interval insertion.");

    UpdateH(insResult.first);
    UpdateH(--insResult.first);

    if(!mNeedRefillQueue)
    {
      pNewInterval->R = CalculateR(pNewInterval);
      mNextIntervals[i]->R = CalculateR(mNextIntervals[i]);
      mQueue.push(pNewInterval);
      mQueue.push(pOldInterval);
    }
  }
}

void NLPSolver::CalculateNextPoints()
{
  for(size_t i = 0; i < mParameters.numThreads; i++)
  {
    mNextIntervals[i] = mQueue.top();
    mQueue.pop();
    mNextPoints[i].x = GetNextPointCoordinate(mNextIntervals[i]);

    if (mNextPoints[i].x >= mNextIntervals[i]->pr.x || mNextPoints[i].x <= mNextIntervals[i]->pl.x)
      throw std::runtime_error("The next point is outside of the subdivided interval");

    mEvolvent.GetImage(mNextPoints[i].x, mNextPoints[i].y);
  }
}

void NLPSolver::RefillQueue()
{
  mQueue = PriorityQueue();
  for (const auto& pInterval : mSearchInformation)
  {
    pInterval->R = CalculateR(pInterval);
    mQueue.push(pInterval);
  }
  mNeedRefillQueue = false;
}

void NLPSolver::EstimateOptimum()
{
  for (size_t i = 0; i < mNextPoints.size(); i++)
  {
    if (mOptimumEstimation.idx < mNextPoints[i].idx ||
        mOptimumEstimation.idx == mNextPoints[i].idx &&
        mOptimumEstimation.g[mOptimumEstimation.idx] > mNextPoints[i].g[mNextPoints[i].idx])
    {
      mOptimumEstimation = mNextPoints[i];
    }
  }
}

void NLPSolver::UpdateH(std::set<Interval*, CompareIntervals>::iterator it_iter)
{
  /*
  Trial& currentPoint = mSearchData[idx];
  int left_idx = idx - 1;
  while(left_idx > 0 && mSearchData[left_idx].v != currentPoint.v)
    left_idx--;
  if(left_idx != (int)idx && mSearchData[left_idx].v == mSearchData[idx].v)
    UpdateMu(mSearchData[left_idx], mSearchData[idx]);

  size_t right_idx = idx + 1;
  while(right_idx < mSearchData.size() - 1 && mSearchData[right_idx].v != currentPoint.v)
    right_idx++;
  if(right_idx != idx && mSearchData[right_idx].v == mSearchData[idx].v)
    UpdateMu(mSearchData[idx], mSearchData[right_idx]);
  }
  */
  Interval* it = *it_iter;
  if (it->pr.idx != it->pl.idx)
    return;
  double oldMu = mHEstimations[it->pl.idx];
  double newMu = fabs(it->pr.g[it->pr.idx] - it->pl.g[it->pl.idx]) /
    pow(it->pr.x - it->pl.x, 1. / mProblem->GetDimension());

  if (newMu > oldMu || (oldMu == 1.0 && newMu > zeroHLevel))
  {
    mHEstimations[it->pr.idx] = newMu;
    mNeedRefillQueue = true;
  }
}

double NLPSolver::CalculateR(Interval* i) const
{
  if(i->pl.idx == i->pr.idx)
  {
    const int v = i->pr.idx;
    return i->delta + pow((i->pr.g[v] - i->pl.g[v]) / (mParameters.r * mHEstimations[v]), 2) / i->delta -
      2.*(i->pr.g[v] + i->pl.g[v] - 2*mZEstimations[v]) / (mParameters.r * mHEstimations[v]);
  }
  else if(i->pl.idx < i->pr.idx)
    return 2*i->delta - 4*(i->pr.g[i->pr.idx] - mZEstimations[i->pr.idx]) / (mParameters.r * mHEstimations[i->pr.idx]);
  else
    return 2*i->delta - 4*(i->pl.g[i->pl.idx] - mZEstimations[i->pl.idx]) / (mParameters.r * mHEstimations[i->pl.idx]);
}

double NLPSolver::GetNextPointCoordinate(Interval* i) const
{
  double x;
  if(i->pr.idx == i->pl.idx)
  {
    const int v = i->pr.idx;
    double dg = i->pr.g[v] - i->pl.g[v];
    x = 0.5 * (i->pr.x + i->pl.x) -
      0.5*((dg > 0.) ? 1. : -1.) * pow(fabs(dg) / mHEstimations[v], mProblem->GetDimension()) / mParameters.r;
  }
  else
    x = 0.5 * (i->pr.x + i->pl.x);

  if (x >= i->pr.x || x <= i->pl.x)
    throw std::runtime_error("Point is outside of interval\n");

  return x;
}
