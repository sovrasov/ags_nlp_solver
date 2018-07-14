#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <memory>

#include <gkls_function.hpp>
#include <grishagin_function.hpp>
#include <cmdline.h>

#include "solver.hpp"

using namespace ags;

void saveStatistics(const std::vector<std::vector<unsigned>>& stat, const cmdline::parser& parser);
void initParser(cmdline::parser& parser);

int main(int argc, char** argv)
{
  cmdline::parser parser;
  initParser(parser);
  parser.parse_check(argc, argv);

  SolverParameters parameters;
  parameters.eps = parser.get<double>("accuracy");
  parameters.r = parser.get<double>("reliability");
  parameters.itersLimit = parser.get<int>("itersLimit");
  parameters.evolventDensity = parser.get<int>("evolventDensity");
  parameters.refineSolution = parser.exist("refineLoc");

  std::string problemClass = parser.get<std::string>("problemsClass");

  auto start = std::chrono::system_clock::now();
  std::vector<std::vector<unsigned>> allStatistics;

  double objectiveAvgConst = 0.;
  double solutionCheckAcc = 3*parser.get<double>("accuracy");

//#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < 100; i++)
  {
    std::shared_ptr<IGOProblem<double>> problem;
    NLPSolver solver;
    solver.SetParameters(parameters);
    if (problemClass == "gklsS" || problemClass == "gklsH")
    {
      auto *func = new gkls::GKLSFunction();
      if (problemClass == "gklsS")
        func->SetFunctionClass(gkls::Simple, parser.get<int>("dim"));
      else
        func->SetFunctionClass(gkls::Hard, parser.get<int>("dim"));

      func->SetType(gkls::TD);
      func->SetFunctionNumber(i + 1);
      problem = std::shared_ptr<IGOProblem<double>>(func);
    }
    else if (problemClass == "grish")
    {
      auto *func = new vagrish::GrishaginFunction();
      func->SetFunctionNumber(i + 1);
      problem = std::shared_ptr<IGOProblem<double>>(func);
    }

    solver.SetProblem(problem);
    Trial optimalPoint;
    try
    {
      optimalPoint = solver.Solve();
    }
    catch (const std::runtime_error& err)
    {
      std::cout << "Exception in solver! " << std::string(err.what()) << "\n";
    }
//#pragma omp critical
    {
      allStatistics.push_back(solver.GetCalculationsStatistics());
      objectiveAvgConst += solver.GetHolderConstantsEstimations().back();

      double optPoint[solverMaxDim];
      problem->GetOptimumPoint(optPoint);
      bool isSolved = !solver_utils::checkVectorsDiff(
        optPoint, optimalPoint.y, problem->GetDimension(), solutionCheckAcc);
      std::cout << "Problem # " << i + 1;
      if (isSolved)
      {
        std::cout << " solved.";
        allStatistics.back().push_back(1);
      }
      else
      {
        std::cout << " not solved.";
        allStatistics.back().push_back(0);
      }
      std::cout << " Iterations performed: " << allStatistics.back()[0] << "\n";
    }
  }
  auto end = std::chrono::system_clock::now();
  objectiveAvgConst /= 100;

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Time elapsed: " << elapsed_seconds.count() << "s\n";
  std::cout << "Objective average Holder const estimation: " << objectiveAvgConst << "\n";

  saveStatistics(allStatistics, parser);

  return 0;
}

void saveStatistics(const std::vector<std::vector<unsigned>>& stat, const cmdline::parser& parser)
{
  size_t numFuncs = stat.back().size() - 1;
  std::vector<double> avgCalcs(numFuncs, 0.);
  unsigned solvedCounter = 0;
  unsigned maxIters = 0;

  for(const auto& elem : stat)
  {
    maxIters = std::max(maxIters, elem[0]);
    for(size_t j = 0; j < numFuncs; j++)
      avgCalcs[j] += elem[j];
    solvedCounter += elem.back();
  }
  for(size_t j = 0; j < numFuncs; j++)
  {
    avgCalcs[j] /= stat.size();
    std::cout << "Average calculations number of function # " << j << " = " << avgCalcs[j] << "\n";
  }
  std::cout << "Problems solved: " << solvedCounter << "\n";
  std::cout << "Maximum number of iterations: " << maxIters << "\n";

  if(parser.exist("saveStat"))
  {
    std::vector<std::pair<int, int>> operationCharacteristic;
    const unsigned opStep = maxIters / 150;
    for(unsigned i = 0; i < maxIters + opStep; i+= opStep)
    {
      int solvedProblemsCnt = 0;
      for(const auto& elem : stat)
        if(elem.back() && elem[numFuncs - 1] <= i)
          solvedProblemsCnt++;
      operationCharacteristic.push_back(std::make_pair(i, solvedProblemsCnt));
    }

    auto fileName = parser.get<std::string>("outFile");
    const std::string sep = "_";
    const std::string stopType = parser.exist("accuracyStop") ? "accuracy" : "optPoint";
    std::string generatedName = parser.get<std::string>("problemsClass") + sep +
      "n_" + std::to_string(parser.get<int>("dim")) + sep +
      "r_" + std::to_string(parser.get<double>("reliability")) + sep +
      "eps_" + std::to_string(parser.get<double>("accuracy")) + sep +
      "lm_" + std::to_string(parser.get<int>("localMix"));
    if(fileName.empty())
      fileName = generatedName + ".csv";

    std::ofstream fout;
    fout.open(fileName, std::ios_base::out);
    fout << generatedName << std::endl;
    for(const auto& point : operationCharacteristic)
      fout << point.first << ", " << point.second << std::endl;
  }
}

void initParser(cmdline::parser& parser)
{
  parser.add<int>("evolventDensity", 'm', "", false, 12,
    cmdline::range(9, 16));
  parser.add<double>("reliability", 'r', "reliability parameter for the method",
    false, 5, cmdline::range(1., 1000.));
  parser.add<double>("accuracy", 'e', "accuracy of the method", false, 0.01);
  parser.add<double>("reserves", 'E', "eps-reserves for all constraints", false, 0);
  parser.add<int>("itersLimit", 'i', "limit of iterations for the method", false, 5000);
  parser.add<int>("dim", 'd', "test problem dimension (will be set if supported)", false, 2);
  parser.add<int>("localMix", 'q', "local mix parameter", false, 0, cmdline::range(-20, 20));
  parser.add<std::string>("problemsClass", 'c', "Name of the used problems class", false,
    "gklsS", cmdline::oneof<std::string>("gklsS", "gklsH", "grish"));
  parser.add<std::string>("outFile", 'f', "Name of the output .csv file with statistics", false,
    "");
  parser.add("saveStat", 's', "Save statistics in a .csv file");
  parser.add("refineLoc", 'l', "Refine the global solution using a local optimizer");
}
