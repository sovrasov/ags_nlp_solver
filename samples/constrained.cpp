#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

#include "solver.hpp"

int main(int argc, char* argv[])
{
  auto parameters = SolverParameters(0.001, 3., 1, 10000);
  parameters.rEps = 0.01;

  NLPSolver solver;
  solver.SetParameters(parameters);
  solver.SetProblem({
    [](const double* x) {return 0.01*(pow(x[0] - 2.2, 2) + pow(x[1] - 1.2, 2) - 2.25);},
    [](const double* x) {return 100 * (1 - pow(x[0] - 2, 2) / 1.44 - pow(0.5*x[1], 2));},
    [](const double* x) {return 10 * (x[1] - 1.5 - 1.5*sin(2*M_PI*(x[0] - 1.75)));},
    [](const double* x) {return -1.5*pow(x[0], 2) * exp(1 - pow(x[0], 2)
				- 20.25*pow(x[0] - x[1], 2)) - pow(0.5 * (x[1] - 1)*(x[0]- 1), 4)
				* exp(2 - pow(0.5 * (x[0] - 1), 4) - pow(x[1] - 1, 4));}
  }, {0, -1}, {4, 3});

  auto optimalPoint = solver.Solve();
  auto calcCounters = solver.GetCalculationsStatistics();
  auto holderConstEstimations = solver.GetHolderConstantsEstimations();

  for (size_t i = 0; i < calcCounters.size() - 1; i++)
    std::cout << "Number of calculations of constraint # " << i << ": " << calcCounters[i] << "\n";
  std::cout << "Number of calculations of objective: " << calcCounters.back() << "\n";

  for (size_t i = 0; i < holderConstEstimations.size() - 1; i++)
    std::cout << "Estimation of Holder constant of function # " << i << ": " << holderConstEstimations[i] << "\n";
  std::cout << "Estimation of Holder constant of objective: " << holderConstEstimations.back() << "\n";

  std::cout << "Optimal value: " << optimalPoint.g[optimalPoint.idx] << "\n";
  std::cout << "x = " << optimalPoint.y[0] << " y = " << optimalPoint.y[1] << "\n";

  return 0;
}
