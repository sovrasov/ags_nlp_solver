#include "solver.hpp"

NLPSolver::NLPSolver() {}

void NLPSolver::SetParameters(const SolverParameters& params)
{

}

void NLPSolver::SetProblem(std::shared_ptr<IGOProblem<double>> problem)
{

}

Trial NLPSolver::Solve()
{
    return Trial();
}

std::vector<int> NLPSolver::GetCalculationsStatistics() const
{
    return std::vector<int>(1, 0);
}

std::vector<double> NLPSolver::GetHolderConstantsEstimations() const
{
    return std::vector<double>(1, 0);
}
