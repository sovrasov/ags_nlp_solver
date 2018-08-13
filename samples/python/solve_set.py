import ags_solver
import go_problems
from benchmark_tools.core import Solver, solve_class, GrishClass, GKLSClass
from benchmark_tools.plot import plot_cmcs
from benchmark_tools.stats import save_stats, compute_stats

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
    params.eps = 0.001
    solver = ags_solver.Solver()
    solver.SetParameters(params)

    problems = GrishClass()
    calc_stats, solved_status = solve_class(problems, AGSWrapper(solver))
    cmc_curve, avg_calculations, num_solved = compute_stats(calc_stats, solved_status)

    print('Problems solved: {}'.format(num_solved))
    for i, avg in enumerate(avg_calculations[:-1]):
        print('Average number of calculations of constraint #{}: {}'.format(i, avg))
    print('Average number of calculations of objective: {}'.format(avg_calculations[-1]))

    plot_cmcs([cmc_curve], show=True, filename='')
    #save_stats(cmc_curve, avg_calculations)

if __name__ == '__main__':
    main()
