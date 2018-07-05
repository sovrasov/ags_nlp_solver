#include <iostream>

#include "solver.hpp"

int main(int argc, char* argv[])
{
  auto parameters = SolverParameters(0.01, 3., 1, 5000);


  NLPSolver solver;
  solver.SetParameters(parameters);
  solver.SetProblem({
    [](const double* y) {return (y[0]-.5)*(y[0]-.5) + y[1]*y[1] - 0.15;},
    [](const double* y) {return y[0]*y[0] + y[1]*y[1];}
  }, {-1, -1}, {1, 1});

  auto optimalPoint = solver.Solve();
  auto calcCounters = solver.GetCalculationsStatistics();

  for (size_t i = 0; i < calcCounters.size() - 1; i++)
    std::cout << "Number of calculations of constraint # " << i << ": " << calcCounters[i] << "\n";
  std::cout << "Number of calculations of objective: " << calcCounters.back() << "\n";

  std::cout << "Optimal value: " << optimalPoint.g[optimalPoint.idx] << "\n";
  std::cout << "x = " << optimalPoint.y[0] << " y = " << optimalPoint.y[1] << "\n";

  return 0;
}
