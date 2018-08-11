import ags_solver
import math
import argparse

def c1(x):
    return 0.01*(math.pow(x[0] - 2.2, 2) + math.pow(x[1] - 1.2, 2) - 2.25)
def c2(x):
    return 100 * (1 - math.pow(x[0] - 2, 2) / 1.44 - math.pow(0.5*x[1], 2))
def c3(x):
    return 10 * (x[1] - 1.5 - 1.5*math.sin(2*math.pi*(x[0] - 1.75)))
def obj_f(x):
    return -1.5*math.pow(x[0], 2) * math.exp(1 - pow(x[0], 2)
    - 20.25*math.pow(x[0] - x[1], 2)) - math.pow(0.5 * (x[1] - 1)*(x[0]- 1), 4) \
    * math.exp(2 - math.pow(0.5 * (x[0] - 1), 4) - math.pow(x[1] - 1, 4))

bounds = [[0, -1], [4, 3]]

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
    solver.SetProblem([c1, c2, c3, obj_f], *bounds)
    point, val, idx = solver.Solve()

    calcCounters = solver.GetCalculationsStatistics()
    holderConstEstimations = solver.GetHolderConstantsEstimations()

    for i in range(len(calcCounters) - 1):
        print('Number of calculations of constraint # {}: {}'.format(i, calcCounters[i]))
    print('Number of calculations of objective: {}'.format(calcCounters[-1]))

    for i in range(len(holderConstEstimations) - 1):
        print('Estimation of Holder constant of function # {}: {}'.format(i, holderConstEstimations[i]))
    print('Estimation of Holder constant of objective: {}'.format(holderConstEstimations[-1]))

    if idx < 3:
        print('Feasible point not found')
    else:
        print('Optimal value: {}'.format(val))
        print('x = ' + str(point[0]) + ' y = ' + str(point[1]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sample for AGS solver')
    parser.add_argument('--m', type=int, default=12, help='Evolvent density')
    parser.add_argument('--r', type=float, default=3, help='Reliability parameter for the solver')
    parser.add_argument('--eps', type=float, default=0.001, help='Accuracy of the method')
    parser.add_argument('--epsR', type=float, default=0.01, help='eps-reserves for all constraints')
    parser.add_argument('--max_iters', type=int, default=10000, help='limit of iterations for the method')
    parser.add_argument('--refine_loc', action='store_true', help='Refine the global solution using a local optimizer')

    main(parser.parse_args())
