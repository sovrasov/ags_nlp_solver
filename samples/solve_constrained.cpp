#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cmdline.h>

#include "solver.hpp"

using namespace ags;

void initParser(cmdline::parser& parser);

int main(int argc, char* argv[])
{
  cmdline::parser parser;
  initParser(parser);
  parser.parse_check(argc, argv);

  SolverParameters parameters;
  parameters.eps = parser.get<double>("accuracy");
  parameters.r = parser.get<double>("reliability");
  parameters.itersLimit = parser.get<int>("itersLimit");
  parameters.epsR = parser.get<double>("reserves");
  parameters.evolventDensity = parser.get<int>("evolventDensity");
  parameters.refineSolution = parser.exist("refineLoc");
  parameters.localMix = parser.get<int>("localMix");

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

  if (optimalPoint.idx < 3)
    std::cout << "Feasible point not found" << "\n";
  else
  {
    std::cout << "Optimal value: " << optimalPoint.g[optimalPoint.idx] << "\n";
    std::cout << "x = " << optimalPoint.y[0] << " y = " << optimalPoint.y[1] << "\n";
  }

  return 0;
}

void initParser(cmdline::parser& parser)
{
  parser.add<int>("evolventDensity", 'm', "", false, 12,
    cmdline::range(9, 16));
  parser.add<double>("reliability", 'r', "reliability parameter for the method",
    false, 3, cmdline::range(1., 1000.));
  parser.add<double>("accuracy", 'e', "accuracy of the method", false, 0.001);
  parser.add<double>("reserves", 'E', "eps-reserves for all constraints", false, 0.01);
  parser.add<int>("itersLimit", 'i', "limit of iterations for the method", false, 10000);
  parser.add("refineLoc", 'l', "Refine the global solution using a local optimizer");
  parser.add<int>("localMix", 'q', "local mix parameter", false, 0, cmdline::range(-20, 20));
}
