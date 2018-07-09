#pragma once

#include "data_types.hpp"
#include "problem_interface.hpp"

#include <memory>
#include <vector>

class HookeJeevesOptimizer
{
private:
  double mEps;
  double mStep;
  double mStepMultiplier;

  std::shared_ptr<IGOProblem<double>> mProblem;

  Trial mCurrentPoint;
  Trial mStartPoint;
  Trial mCurrentResearchDirection;
  Trial mPreviousResearchDirection;

  void DoStep();
  double ComputeObjective(const double* x) const;
  double MakeResearch(double*);

public:
  void SetParameters(double eps, double step, double stepMult);
	Trial Optimize(std::shared_ptr<IGOProblem<double>> problem, const Trial& startPoint);
};
