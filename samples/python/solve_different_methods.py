import functools
import numpy as np
import argparse
import ags_solver
import go_problems
import nlopt
import sys
sys.path.append('./')
from Simple import SimpleTuner
import itertools
from scipy.spatial import Delaunay

from scipy.optimize import differential_evolution

from benchmark_tools.core import Solver, solve_class, GrishClass, GKLSClass
from benchmark_tools.plot import plot_cmcs
from benchmark_tools.stats import save_stats, compute_stats

class AGSWrapper(Solver):
    def __init__(self, dist_stop, eps=0.01):
        params = ags_solver.Parameters()
        if dist_stop:
            params.eps = 0
            params.r = 4.4
            params.itersLimit = 40000
        self.solver = ags_solver.Solver()
        self.solver.SetParameters(params)
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

class SCDEWrapper(Solver):
    def __init__(self, dist_stop, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps

    def Solve(self, problem):
        if self.dist_stop:
            problem.SetEps(self.eps)
        lb, ub = problem.GetBounds()
        bounds = [(l, u) for l, u in zip(lb, ub)]
        result = \
            differential_evolution(
            lambda x: problem.Calculate(x), bounds, mutation=(1.1,1.9),
            tol=1e-3, maxiter=3000, popsize=60, disp=False, seed=100)

        n_evals = problem.GetCalculationsStatistics()
        return result.x, result.fun, n_evals

class PyEvolveWrapper(Solver):
    def __init__(self, dist_stop, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps

    def Solve(self, problem):
        if self.dist_stop:
            problem.SetEps(self.eps)
        lb, ub = problem.GetBounds()

        # Genome instance
        genome = G1DList.G1DList(2)
        genome.setParams(rangemin=lb[0], rangemax=ub[0], bestRawScore=-100, roundDecimal=2)
        genome.initializator.set(Initializators.G1DListInitializatorReal)
        genome.mutator.set(Mutators.G1DListMutatorRealGaussian)

        # The evaluator function (objective function)
        genome.evaluator.set(lambda x: problem.Calculate(x) + 100)

        # Genetic Algorithm Instance
        ga = GSimpleGA.GSimpleGA(genome)
        ga.selector.set(Selectors.GRouletteWheel)

        ga.minimax = Consts.minimaxType["minimize"]
        ga.setGenerations(5000)
        ga.setMutationRate(0.05)
        ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

        # Do the evolution, with stats dump
        # frequency of 10 generations
        ga.evolve(freq_stats=100)

        # Best individual
        best = ga.bestIndividual()
        print ("\nBest individual score: %.2f" % (best.score - 100,))
        print (best)

from bayes_opt import BayesianOptimization
class BOptWrapper:
    def __init__(self, dist_stop, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps

    def Solve(self, problem):
        lb, ub = problem.GetBounds()
        bo = BayesianOptimization(lambda x, y: -problem.Calculate([x, y]),
                          {'x': (lb[0], ub[0]), 'y': (lb[1], ub[1])})
        bo.maximize(init_points=5, n_iter=20, kappa=1.5)
        n_evals = problem.GetCalculationsStatistics()
        opt_val = -bo.res['max']['max_val']
        opt_point = [bo.res['max']['max_params']['x'], bo.res['max']['max_params']['y']]
        return opt_point, opt_val, n_evals

class SimpleWrapper:
    def __init__(self, dist_stop, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps

    def Solve(self, problem):
        if self.dist_stop:
            problem.SetEps(self.eps)
        objective_function = lambda x: -problem.Calculate(x)
        lb, ub = problem.GetBounds()
        bounds = [[l, u] for l, u in zip(lb, ub)]
        points = np.array([point for point in itertools.product(*bounds)])
        tri = Delaunay(points)
        optimization_domain_vertices = points[tri.simplices]
        number_of_iterations = 40000
        exploration = 0.1 # optional, default 0.15
        tuner = SimpleTuner(optimization_domain_vertices, objective_function, exploration_preference=exploration)
        tuner.optimize(number_of_iterations)
        opt_val, opt_point = tuner.get_best()
        print(-opt_val)
        #tuner.plot() # only works in 2D
        n_evals = problem.GetCalculationsStatistics()
        return opt_point, -opt_val, n_evals

class NLOptWrapper:
    def __init__(self, dist_stop, method=nlopt.GD_STOGO, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.method = method

    def Solve(self, problem):
        if self.dist_stop:
            problem.SetEps(self.eps)
        lb, ub = problem.GetBounds()

        self.opt = nlopt.opt(self.method, problem.GetDimension())
        self.opt.set_local_optimizer(nlopt.opt(nlopt.LN_SBPLX, problem.GetDimension()))
        self.opt.set_lower_bounds(lb)
        self.opt.set_upper_bounds(ub)
        self.opt.set_min_objective(lambda x, grad: problem.Calculate(x))
        self.opt.set_maxeval(60000)
        self.opt.set_xtol_rel(1e-13)
        #self.opt.set_population(5000)

        x = self.opt.optimize([.5]*problem.GetDimension())
        minf = self.opt.last_optimum_value()

        n_evals = problem.GetCalculationsStatistics()
        return x, minf, n_evals

algos = {'scd': SCDEWrapper, 'ags': AGSWrapper,
         'direct': functools.partial(NLOptWrapper, method=nlopt.GN_ORIG_DIRECT),
         'directl': functools.partial(NLOptWrapper, method=nlopt.GN_ORIG_DIRECT_L),
         'stogo': functools.partial(NLOptWrapper, method=nlopt.GD_STOGO),
         'mlsl': functools.partial(NLOptWrapper, method=nlopt.G_MLSL_LDS),
         'crs': functools.partial(NLOptWrapper, method=nlopt.GN_CRS2_LM),
         'simple': SimpleWrapper}

def algo2cature(algo):
    if algo == 'scd':
        return 'Scipy DE'
    if algo == 'ags':
        return 'AGS'
    if algo == 'direct':
        return 'DIRECT'
    if algo == 'directl':
        return 'DIRECTl'
    if algo == 'simple':
        return 'Simple'
    if algo == 'stogo':
        return 'StoGO'
    if algo == 'mlsl':
        return 'MLSL'
    if algo == 'crs':
        return 'CRS'

def main(args):

    wrapper = algos[args.algo]

    if args.problems_class == 'grish':
        problems = GrishClass()
    else:
        assert args.problems_dim > 1 and args.problems_dim < 6
        if args.problems_class == 'gklss':
            problems = GKLSClass(args.problems_dim, go_problems.GKLSClass.Simple)
        else:
            problems = GKLSClass(args.problems_dim, go_problems.GKLSClass.Hard)

    calc_stats, solved_status = solve_class(problems, wrapper(args.dist_stop), verbose=args.verbose, eps_check=0.01)
    stats = compute_stats(calc_stats, solved_status)

    print('Problems solved: {}'.format(stats['num_solved']))
    for i, avg in enumerate(stats['avg_calcs'][:-1]):
        print('Average number of calculations of constraint #{}: {}'.format(i, avg))
    print('Average number of calculations of objective: {}'.format(stats['avg_calcs'][-1]))

    plot_cmcs([stats['cmc']], captures=[algo2cature(args.algo)], show=True, filename='')
    save_stats(stats, args.stats_fname, capture=algo2cature(args.algo))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample for AGS solver')
    parser.add_argument('--max_iters', type=int, default=10000, help='limit of iterations for the method')
    parser.add_argument('--problems_class', type=str, choices=['grish','gklss','gklsh'], default='grish')
    parser.add_argument('--algo', type=str, choices=algos.keys(), default='scd')
    parser.add_argument('--problems_dim', type=int, default=2)
    parser.add_argument('--verbose', action='store_true', help='Print additional info to console')
    parser.add_argument('--dist_stop', action='store_true', help='Stop algorithm then the next point is close enough to the optimum')
    parser.add_argument('--stats_fname', type=str, default='')

    main(parser.parse_args())
