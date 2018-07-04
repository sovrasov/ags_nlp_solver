#pragma once

#include "data_types.hpp"
#include "evolvent.hpp"
#include "problem_interface.hpp"

#include <vector>
#include <memory>
#include <queue>
#include <set>

struct SolverParameters
{
  double eps;
  double r;
  unsigned numThreads;
  unsigned trialsLimit;
  unsigned evolventTightness = 12;
  double rEps = 0.;

  SolverParameters() {}
  SolverParameters(double _eps, double _r,
      unsigned _numThreads, unsigned _trialsLimit) :
        eps(_eps), r(_r), numThreads(_numThreads), trialsLimit(_trialsLimit)
  {}
};

class NLPSolver
{
protected:
  using PriorityQueue =
    std::priority_queue<Interval*, std::vector<Interval*>, CompareByR>;

  SolverParameters mParameters;
  std::shared_ptr<IGOProblem<double>> mProblem;
  Evolvent mEvolvent;

  std::vector<double> mHEstimations;
  std::vector<double> mZEstimations;
  std::vector<Trial> mNextPoints;
  PriorityQueue mQueue;
  std::set<Interval*> mSearchInformation;
  std::vector<Interval*> mNextIntervals;
  Trial mOptimumEstimation;

  std::vector<unsigned> mCalculationsCounters;
  unsigned mIterationsCounter;
  bool mNeedRefillQueue;
  double mMinDelta;
  int mMaxIdx;

  void FirstIteration();
  void MakeTrials();
  void InsertIntervals();
  void CalculateNextPoints();
  void RefillQueue();
  void EstimateOptimum();

  void InitDataStructures();
  void ClearDataStructures();

  void UpdateH(std::set<Interval*>::iterator) const;
  double CalculateR(Interval*) const;
  double GetNextPointCoordinate(Interval*) const;

public:
  NLPSolver();

  void SetParameters(const SolverParameters& params);
  void SetProblem(std::shared_ptr<IGOProblem<double>> problem);

  Trial Solve();
  std::vector<unsigned> GetCalculationsStatistics() const;
  std::vector<double> GetHolderConstantsEstimations() const;
};

namespace solver_utils
{
  inline bool checkVectorsDiff(const double* y1, const double* y2, size_t dim, double eps)
  {
    for (size_t i = 0; i < dim; i++)
    {
      if (fabs(y1[i] - y2[i]) > eps)
        return true;
    }

    return false;
  }
}
