import functools
import numpy as np
import math
import argparse
import ags_solver
import go_problems
import nlopt
import sys
from Simple import SimpleTuner
import itertools
from scipy.spatial import Delaunay
from scipy.optimize import differential_evolution
from scipy.optimize import basinhopping
from sdaopt import sda
from stochopy import Evolutionary
from pySOT import *
from poap.controller import SerialController, BasicWorkerThread, ThreadController
from pyOpt import Optimization
from pyOpt import MIDACO
import pyOpt
import pygmo

from benchmark_tools.core import Solver, solve_class, GrishClass, GKLSClass
from benchmark_tools.plot import plot_cmcs
from benchmark_tools.stats import save_stats, compute_stats

class AGSWrapper(Solver):
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        params = self.class_name2params(class_name)
        if dist_stop:
            params.eps = 0
            params.itersLimit = max_iters
        self.solver = ags_solver.Solver()
        self.solver.SetParameters(params)
        self.dist_stop = dist_stop
        self.eps = eps

    def class_name2params(self, name):
        params = ags_solver.Parameters()
        if 'grish' in name:
            params.r = 3
        elif 'gklss2' in name:
            params.r = 4.6
        elif 'gklsh2' in name:
            params.r = 6.5
        elif 'gklss3' in name:
            params.r = 3.7
        elif 'gklsh3' in name:
            params.r = 4.4
        elif 'gklss4' in name:
            params.r = 4.7
        elif 'gklsh4' in name:
            params.r = 4.9
        elif 'gklss5' in name:
            params.r = 4
            params.evolventDensity = 10
        elif 'gklsh5' in name:
            params.r = 4
            params.evolventDensity = 10
        return params

    def Solve(self, problem):
        self.solver.SetProblem([lambda x: problem.Calculate(x)], *problem.GetBounds())
        #self.solver.SetProblem(problem)
        if not self.dist_stop:
            point, val, idx = self.solver.Solve()
        else:
            opt_pt = np.array(problem.GetOptimumPoint())
            point, val, idx = self.solver.Solve(lambda x: np.linalg.norm(np.array(x)-opt_pt, np.inf) < self.eps)
        #calcCounters = self.solver.GetCalculationsStatistics()
        calcCounters = problem.GetCalculationsStatistics()
        return point, val, calcCounters

class SDAWrapper:
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.max_iters = max_iters
        self.class_name = class_name
    def Solve(self, problem):
        lb, ub = problem.GetBounds()
        ret = sda(lambda x: problem.Calculate(x), None, bounds=list(zip(lb, ub)), \
                  seed=100, maxfun=self.max_iters, visit=2.72, maxiter=self.max_iters)
        n_evals = problem.GetCalculationsStatistics()
        return ret.x, ret.fun, n_evals

class SCBasinhoppingWrapper:
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.max_iters = max_iters
        self.class_name = class_name

    def Solve(self, problem):
        lb, ub = problem.GetBounds()
        #pop_size = self.class_name2params(self.class_name)
        class MyBounds(object):
            def __init__(self, xmax=[1.1,1.1], xmin=[-1.1,-1.1] ):
                self.xmax = np.array(xmax)
                self.xmin = np.array(xmin)
            def __call__(self, **kwargs):
                x = kwargs["x_new"]
                tmax = bool(np.all(x <= self.xmax))
                tmin = bool(np.all(x >= self.xmin))
                return tmax and tmin

        x0 = [.5]*problem.GetDimension()
        result = \
            basinhopping(lambda x: problem.Calculate(x), x0, accept_test=MyBounds(ub, lb), seed=100, T=10, stepsize=0.3)

        n_evals = problem.GetCalculationsStatistics()
        return result.x, result.fun, n_evals

class SCDEWrapper(Solver):
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.max_iters = max_iters
        self.class_name = class_name

    def class_name2params(self, name):
        if 'grish' in name:
            popsize = 60
        elif 'gklss2' in name:
            popsize = 60
        elif 'gklsh2' in name:
            popsize = 60
        elif 'gklss3' in name:
            popsize = 70
        elif 'gklsh3' in name:
            popsize = 80
        elif 'gklss4' in name:
            popsize = 90
        elif 'gklsh4' in name:
            popsize = 100
        elif 'gklss5' in name:
            popsize = 120
        elif 'gklsh5' in name:
            popsize = 140
        return popsize

    def Solve(self, problem):
        lb, ub = problem.GetBounds()
        bounds = [(l, u) for l, u in zip(lb, ub)]
        pop_size = self.class_name2params(self.class_name)
        result = \
            differential_evolution(
            lambda x: problem.Calculate(x), bounds, mutation=(1.1,1.9),
            tol=1e-12, maxiter=int(float(self.max_iters) / (pop_size*problem.GetDimension())), popsize=pop_size, disp=False, seed=100)

        n_evals = problem.GetCalculationsStatistics()
        return result.x, result.fun, n_evals

class PyEvolveWrapper(Solver):
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps

    def Solve(self, problem):
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
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
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
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.max_iters = max_iters
        self.exploration = self.class_name2params(class_name)

    def class_name2params(self, name):
        if 'grish' in name:
            return 0.1
        elif 'gklss2' in name:
            return 0.15
        elif 'gklsh2' in name:
            return 0.15
        elif 'gklss3' in name:
            return 0.15
        elif 'gklsh3' in name:
            return 0.25
        elif 'gklss4' in name:
            return 0.2
        elif 'gklsh4' in name:
            return 0.25

    def Solve(self, problem):
        objective_function = lambda x: -problem.Calculate(x)
        lb, ub = problem.GetBounds()
        opt_pt = problem.GetOptimumPoint()
        bounds = [[l, u] for l, u in zip(lb, ub)]
        points = np.array([point for point in itertools.product(*bounds)])
        tri = Delaunay(points)
        optimization_domain_vertices = points[tri.simplices]
        exploration = self.exploration # optional, default 0.15
        tuner = SimpleTuner(optimization_domain_vertices, objective_function, \
                exploration_preference=exploration,
                stop_criterion=lambda x:np.linalg.norm(np.array(x)-opt_pt, np.inf) < self.eps)
        tuner.optimize(self.max_iters)
        opt_val, opt_point = tuner.get_best()
        #tuner.plot() # only works in 2D
        n_evals = problem.GetCalculationsStatistics()
        return opt_point, -opt_val, n_evals

class NLOptWrapper:
    def __init__(self, dist_stop, max_iters, class_name, method=nlopt.GD_STOGO, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.method = method
        self.max_iters = max_iters
        self.pop_size = self.class_name2params(class_name)

    def class_name2params(self, name):
        if 'grish' in name:
            popsize = 150
        elif 'gklss2' in name:
            popsize = 200
        elif 'gklsh2' in name:
            popsize = 400
        elif 'gklss3' in name:
            popsize = 1000
        elif 'gklsh3' in name:
            popsize = 2000
        elif 'gklss4' in name:
            popsize = 8000
        elif 'gklsh4' in name:
            popsize = 16000
        elif 'gklss5' in name:
            popsize = 25000
        elif 'gklsh5' in name:
            popsize = 30000
        return popsize

    def Solve(self, problem):
        lb, ub = problem.GetBounds()
        self.opt = nlopt.opt(self.method, problem.GetDimension())
        self.opt.set_local_optimizer(nlopt.opt(nlopt.LN_SBPLX, problem.GetDimension()))
        self.opt.set_lower_bounds(lb)
        self.opt.set_upper_bounds(ub)
        self.opt.set_min_objective(lambda x, grad: problem.Calculate(x))
        self.opt.set_maxeval(self.max_iters)
        self.opt.set_xtol_rel(1e-13)
        if self.method == nlopt.GN_CRS2_LM:
            self.opt.set_population(self.pop_size)
        x = self.opt.optimize([.5]*problem.GetDimension())
        minf = self.opt.last_optimum_value()
        n_evals = problem.GetCalculationsStatistics()
        return x, minf, n_evals

class StochOpyWrapper:
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.max_iters = max_iters
        self.popsize = self.class_name2params(class_name)

    def class_name2params(self, name):
        if 'grish' in name:
            popsize = 60
        elif 'gklss2' in name:
            popsize = 60
        elif 'gklsh2' in name:
            popsize = 60
        elif 'gklss3' in name:
            popsize = 70
        elif 'gklsh3' in name:
            popsize = 80
        elif 'gklss4' in name:
            popsize = 90
        elif 'gklsh4' in name:
            popsize = 100
        elif 'gklss5' in name:
            popsize = 120
        elif 'gklsh5' in name:
            popsize = 140
        return popsize

    def Solve(self, problem):
        objective_function = lambda x: 50 + problem.Calculate(x)
        lb, ub = problem.GetBounds()
        ea = Evolutionary(objective_function, lower=lb, upper=ub, popsize=self.popsize, \
            max_iter=int(self.max_iters/self.popsize), eps1=1e-16, eps2=1e-16)
        xopt, gfit = ea.optimize(solver='cpso', sync=False, CR=0.4, F=0.5)
        n_evals = problem.GetCalculationsStatistics()
        return xopt, gfit, n_evals

class PySOTWrapper:
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.max_iters = max_iters
        self.popsize = self.class_name2params(class_name)

    def class_name2params(self, name):
        if 'grish' in name:
            popsize = 60
        elif 'gklss2' in name:
            popsize = 60
        elif 'gklsh2' in name:
            popsize = 60
        elif 'gklss3' in name:
            popsize = 70
        elif 'gklsh3' in name:
            popsize = 80
        elif 'gklss4' in name:
            popsize = 90
        elif 'gklsh4' in name:
            popsize = 100
        elif 'gklss5' in name:
            popsize = 120
        elif 'gklsh5' in name:
            popsize = 140
        return popsize

    def Solve(self, problem):
        class pysot_function:
            def __init__(self, problem):
                self.xlow, self.xup = np.array(problem.GetBounds())
                self.dim = problem.GetDimension()
                self.integer = np.array([0])
                self.continuous = np.arange(1, self.dim)
            def objfunction(self, x):
                return problem.Calculate(x)
        objective = pysot_function(problem)
        check_opt_prob(objective)

        maxeval = 1000
        exp_des = SymmetricLatinHypercube(dim=objective.dim, npts=2*objective.dim+1)
        basisp = basis_TD(3, 2)
        surrogate = RBFInterpolant(kernel=CubicKernel, tail=LinearTail, maxp=maxeval)
        adapt_samp = CandidateDYCORS(data=objective, numcand=100*objective.dim)
        controller = SerialController(objective.objfunction)
        #controller = ThreadController()
        nthreads = 1
        strategy = SyncStrategyNoConstraints(
                worker_id=0, data=objective, maxeval=maxeval, nsamples=nthreads,
                exp_design=exp_des, response_surface=surrogate,
                sampling_method=adapt_samp)
        controller.strategy = strategy
        #for _ in range(nthreads):
        #    worker = BasicWorkerThread(controller, objective.objfunction)
        #    controller.launch_worker(worker)
        result = controller.run()
        print(result.value)
        n_evals = problem.GetCalculationsStatistics()
        return result.params[0], result.value, n_evals

class PyOptWrapper:
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.max_iters = max_iters

    def Solve(self, problem):
        objective_function = lambda x: [problem.Calculate(x), 0, 0]
        lb, ub = problem.GetBounds()

        opt_prob = pyOpt.Optimization('Problem', objective_function)
        opt_prob.addObj('f')
        for i in range(problem.GetDimension()):
            opt_prob.addVar('x'+str(i),'c',lower=lb[i],upper=ub[i],value=(lb[i] + ub[i])/2.)
        midaco_none = MIDACO(pll_type=None)
        midaco_none.setOption('IPRINT',-1)
        midaco_none.setOption('ISEED', 100)
        midaco_none.setOption('MAXEVAL',self.max_iters)
        midaco_none.setOption('FOCUS', -4)
        fstr, xstr, inform = midaco_none(opt_prob)

        n_evals = problem.GetCalculationsStatistics()
        return xstr, fstr[0], n_evals
		
 class Pygmo_task:
    func = ''
    bounds = ()
    def __init__(self, lambda_func, lb, ub):
        self.func = lambda_func
        self.bounds = tuple(lb, ub)
        
    def fitness(self, x):
        return [self.func(x)]
    
    def get_bounds(self):
        return self.bounds
		
class PygmoWrapper:
    def __init__(self, dist_stop, max_iters, class_name, eps=0.01):
        self.dist_stop = dist_stop
        self.eps = eps
        self.max_iters = max_iters

    def Solve(self, problem):
        objective_function = lambda x: problem.Calculate(x)
        lb, ub = problem.GetBounds()

        opt_prob = pygmo.problem(Pygmo_task(objective_function, lb, ub))
        
        algo = pygmo.algorithm(pygmo.sade(gen = 1000, variant = 17, variant_adptv = 1, ftol = 0, xtol = 0))
        archi = pygmo.archipelago(3, algo=algo, prob=opt_prob, pop_size=10000)
        archi.evolve()
        archi.wait()

        f = [isl.get_population().champion_f for isl in archi]
        x = [isl.get_population().champion_x for isl in archi]
        import operator
        index_min_f, value_min_f = min(enumerate(f), key=operator.itemgetter(1))
        fstr = str(value_min_f)
        xstr = str(x[index_min_f])
        n_evals = problem.GetCalculationsStatistics()
        
        return xstr, fstr, n_evals

algos = {'scd': SCDEWrapper, 'ags': AGSWrapper,
         'direct': functools.partial(NLOptWrapper, method=nlopt.GN_ORIG_DIRECT),
         'directl': functools.partial(NLOptWrapper, method=nlopt.GN_ORIG_DIRECT_L),
         'stogo': functools.partial(NLOptWrapper, method=nlopt.GD_STOGO),
         'mlsl': functools.partial(NLOptWrapper, method=nlopt.G_MLSL_LDS),
         'crs': functools.partial(NLOptWrapper, method=nlopt.GN_CRS2_LM),
         'simple': SimpleWrapper, 'scb': SCBasinhoppingWrapper,
         'sda': SDAWrapper, 'stochopy': StochOpyWrapper, 'pysot': PySOTWrapper,
         'pyopt': PyOptWrapper,
         'pygmo': PygmoWrapper}

algo2cature = {'scd': 'Scipy DE', 'ags': 'AGS', 'direct': 'DIRECT',
               'directl': 'DIRECTl', 'simple': 'Simple',
               'stogo': 'StoGO', 'mlsl': 'MLSL', 'crs':'CRS', 'scb': 'Scipy B-H',
               'sda': 'SDA', 'stochopy': 'Stochopy', 'pysot': 'PySOT', 'pyopt': 'PyOpt', 'pygmo': 'Pygmo'}

serg_eps = {2: 0.01, 3: 0.01, 4: math.pow(1e-6, 1./4), 5: math.pow(1e-7, 1./5)}

def main(args):

    wrapper_class = algos[args.algo]
    if args.problems_class == 'grish':
        problems = GrishClass()
    else:
        assert args.problems_dim > 1 and args.problems_dim < 6
        if args.problems_class == 'gklss':
            problems = GKLSClass(args.problems_dim, go_problems.GKLSClass.Simple)
        else:
            problems = GKLSClass(args.problems_dim, go_problems.GKLSClass.Hard)

    eps = 0.01
    if args.serg_eps:
        eps = serg_eps[args.problems_dim]
    wrapper = wrapper_class(args.dist_stop, args.max_iters, args.problems_class+str(args.problems_dim), eps=0.01)
    calc_stats, solved_status = solve_class(problems, wrapper, verbose=args.verbose, eps_check=eps)
    stats = compute_stats(calc_stats, solved_status)

    print('Problems solved: {}'.format(stats['num_solved']))
    for i, avg in enumerate(stats['avg_calcs'][:-1]):
        print('Average number of calculations of constraint #{}: {}'.format(i, avg))
    print('Average number of calculations of objective: {}'.format(stats['avg_calcs'][-1]))

    #plot_cmcs([stats['cmc']], captures=[algo2cature(args.algo)], show=True, filename='')
    save_stats(stats, args.stats_fname, capture=algo2cature[args.algo])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample for AGS solver')
    parser.add_argument('--max_iters', type=int, default=10000, help='limit of iterations for the method')
    parser.add_argument('--problems_class', type=str, choices=['grish','gklss','gklsh'], default='grish')
    parser.add_argument('--algo', type=str, choices=algos.keys(), default='scd')
    parser.add_argument('--problems_dim', type=int, default=2)
    parser.add_argument('--verbose', action='store_true', help='Print additional info to console')
    parser.add_argument('--dist_stop', action='store_true', help='Stop algorithm then the next point is close enough to the optimum')
    parser.add_argument('--serg_eps', action='store_true')
    parser.add_argument('--stats_fname', type=str, default='')

    main(parser.parse_args())
