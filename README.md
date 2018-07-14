# Global search NLP solver

An implementation of the algorithm AGS to solve constrained nonlinear programming problems with Lipschitzian functions. AGS was introduced by prof. R.G. Strongin (see R. G. Strongin, D. L. Markin, ,Minimization of multiextremal functions under nonconvex constraints, Cybernetics 22(4), 486-493. Translated from Russian. Consultant Bureau. New York, 1986. [[link]][paper]). The method exploits Peano-type curve to reduce dimension of the source bounded multidimensional constrained NLP problem and then solves a univariate one.

AGS is proven to converge to a global optima if all objectives and constraints satisfy Lipschitz condition in a given hyperrectangle, the reliability parameter `r` is large enough and accuracy parameter `eps` is zero.

## Clone & build, run samples
- on Linux:
```bash
git clone --recursive https://github.com/sovrasov/glob_search_nlp_solver.git
cd multicriterial-go
mkdir build
cd build
cmake ..
make -j 4
./bin/solve_constrained
./bin/solve_set
```
- on Windows:
```batch
git clone --recursive https://github.com/sovrasov/glob_search_nlp_solver.git
cd multicriterial-go
mkdir build
cd build
cmake .. -G "NMake Makefiles"
nmake
.\bin\solve_constrained.exe
.\bin\solve_set.exe
```
[paper]: https://www.tandfonline.com/doi/abs/10.1080/17442508908833568?journalCode=gssr19

## Example of usage
```C++
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

#include "solver.hpp"

using namespace ags;

int main(int argc, char** argv)
{
  auto parameters = SolverParameters();
  parameters.refineSolution = true; // refine solution with a local optimizer

  NLPSolver solver;
  solver.SetParameters(parameters);
  //First 3 functions -- nonlinear inequality constraints g_i(y)<=0
  //Last function -- objective
  //Last 2 arguments -- bounds of the search hyperrectangle
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


  //Optimal point has it's index -- number of the first broken constraint
  //If index equals to the number of constraints, then the point if feasible and
  //objective was evaluated at this point. If the solver returned unfeasible
  //optimal point, the set of feasible points is most likely to be empty.
  if (optimalPoint.idx < 3)
    std::cout << "Feasible point not found" << "\n";
  else
  {
    std::cout << "Optimal value: " << optimalPoint.g[optimalPoint.idx] << "\n";
    std::cout << "x = " << optimalPoint.y[0] << " y = " << optimalPoint.y[1] << "\n";
  }
  return 0;
}
```
