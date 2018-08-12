import ags_solver
import go_problems
from benchmark_tools.core import Solver, solve_class, GrishClass, GKLSClass
from benchmark_tools.plot import plot_cmcs
from benchmark_tools.stats import save_stats, compute_stats

import numpy as np

class AGSWrapper(Solver):
    def __init__(self, solver):
        self.solver = solver

    def Solve(self, problem):
        self.solver.SetProblem(problem)
        point, val, idx = self.solver.Solve()
        calcCounters = self.solver.GetCalculationsStatistics()
        return point, val, calcCounters

def main():

    params = ags_solver.Parameters()
    solver = ags_solver.Solver()
    solver.SetParameters(params)

    problems = GrishClass()
    calc_stats, num_solved = solve_class(problems, AGSWrapper(solver))
    cmc_curve, avg_calculations = compute_stats(calc_stats)
    plot_cmc(cmc_curve)
    save_stats(cmc_curve, avg_calculations)

if __name__ == '__main__':
    main()
