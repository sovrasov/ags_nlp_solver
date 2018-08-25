import numpy as np
import argparse
import ags_solver
import go_problems
from benchmark_tools.core import Solver, solve_class, GrishClass, GKLSClass
from benchmark_tools.plot import plot_cmcs
from benchmark_tools.stats import save_stats, compute_stats

class AGSWrapper(Solver):
    def __init__(self, solver, dist_stop, eps=0.01):
        self.solver = solver
        self.dist_stop = dist_stop
        self.eps = eps

    def Solve(self, problem):
        self.solver.SetProblem(problem)
        if not self.dist_stop:
            point, val, idx = self.solver.Solve()
        else:
            opt_pt = np.array(problem.GetOptimumPoint())
            point, val, idx = self.solver.Solve(lambda x: np.linalg.norm(x-opt_pt, np.inf) < self.eps)
        calcCounters = self.solver.GetCalculationsStatistics()
        return point, val, calcCounters

def main(args):
    params = ags_solver.Parameters()
    params.eps = args.eps
    params.itersLimit = args.max_iters
    params.r = args.r
    params.evolventDensity = args.m
    params.epsR = args.epsR
    params.refineSolution = args.refine_loc

    solver = ags_solver.Solver()
    solver.SetParameters(params)

    if args.problems_class == 'grish':
        problems = GrishClass()
    else:
        assert args.problems_dim > 1 and args.problems_dim < 6
        if args.problems_class == 'gklss':
            problems = GKLSClass(args.problems_dim, go_problems.GKLSClass.Simple)
        else:
            problems = GKLSClass(args.problems_dim, go_problems.GKLSClass.Hard)

    calc_stats, solved_status = solve_class(problems, AGSWrapper(solver, args.dist_stop), verbose=args.verbose)
    stats = compute_stats(calc_stats, solved_status)

    print('Problems solved: {}'.format(stats['num_solved']))
    for i, avg in enumerate(stats['avg_calcs'][:-1]):
        print('Average number of calculations of constraint #{}: {}'.format(i, avg))
    print('Average number of calculations of objective: {}'.format(stats['avg_calcs'][-1]))

    plot_cmcs([stats['cmc']], captures=['AGS'], show=False, filename=args.stats_fname)
    save_stats(stats, args.cmc_fname, capture='AGS')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample for AGS solver')
    parser.add_argument('--m', type=int, default=12, help='Evolvent density')
    parser.add_argument('--r', type=float, default=3, help='Reliability parameter for the solver')
    parser.add_argument('--eps', type=float, default=0.001, help='Accuracy of the method')
    parser.add_argument('--epsR', type=float, default=0.01, help='eps-reserves for all constraints')
    parser.add_argument('--max_iters', type=int, default=10000, help='limit of iterations for the method')
    parser.add_argument('--refine_loc', action='store_true', help='Refine the global solution using a local optimizer')
    parser.add_argument('--stats_fname', type=str, default='AGS_cmc.pdf')
    parser.add_argument('--cmc_fname', type=str, default='ags_cmc.json')
    parser.add_argument('--problems_class', type=str, choices=['grish','gklss','gklsh'], default='grish')
    parser.add_argument('--problems_dim', type=int, default=2)
    parser.add_argument('--verbose', action='store_true', help='Print additional info to console')
    parser.add_argument('--dist_stop', action='store_true', help='Stop algorithm then the next point is close enough to the optimum')

    main(parser.parse_args())
