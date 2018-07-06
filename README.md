# Global search NLP solver

An implementation of the algorithm AGS to solve constrained nonlinear programming problems with Lipschitzian functions. AGS was introduced by prof. R.G. Strongin (see R. G. Strongin, D. L. Markin, ,Minimization of multiextremal functions under nonconvex constraints, Cybernetics 22(4), 486-493. Translated from Russian. Consultant Bureau. New York, 1986. [[link]][paper]). The method exploits Peano-type curve to reduce dimension of the source bounded multidimensional constrained NLP problem and then solves a univariate one.

AGS is proven to converge to a global optima if all objectives and constraints satisfy Lipschitz condition in a given hyperrectangle, the reliability parameter `r` is large enough and accuracy parameter `eps` is zero.

## Clone & build
- on Linux:
```bash
git clone --recursive https://github.com/sovrasov/glob_search_nlp_solver.git
cd multicriterial-go
mkdir build
cd build
cmake ..
make -j 4
```
- on Windows:
```batch
git clone --recursive https://github.com/sovrasov/glob_search_nlp_solver.git
cd multicriterial-go
mkdir build
cd build
cmake .. -G "NMake Makefiles"
nmake
```
[paper]: https://www.tandfonline.com/doi/abs/10.1080/17442508908833568?journalCode=gssr19

## Example of usage
```C++
#include "solver.hpp"

/*.............................*/

auto parameters = SolverParameters(0.001, //tolerance parameter eps
                                   0.01, //reliability parameter r
                                   1, //number of new points per iteration
                                   5000); //max number of iterations
parameters.rEps = 0.01; //

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
```
