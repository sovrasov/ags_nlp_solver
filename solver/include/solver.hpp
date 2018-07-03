#pragma once

#include "data_types.hpp"
#include "problem_interface.hpp"

#include <vector>
#include <memory>

struct SolverParameters
{
  double eps;
  double r;
  unsigned numThreads;
  unsigned trialsLimit;
  unsigned evolventTightness = 12;

  SolverParameters() {}
  SolverParameters(double _eps, double _r,
      unsigned _numThreads, unsigned _trialsLimit) :
        eps(_eps), r(_r), numThreads(_numThreads), trialsLimit(_trialsLimit)
  {}
};


class NLPSolver
{
protected:
public:
  NLPSolver();

  void SetParameters(const SolverParameters& params);
  void SetProblem(std::shared_ptr<IGOProblem<double>> problem);

  Trial Solve();
  std::vector<int> GetCalculationsStatistics() const;
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
