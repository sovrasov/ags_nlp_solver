import ags_solver
import go_problems
from benchmark_tools.core import Solver, solve_class, grish_class, gkls_class
from benchmark_tools.plot import plot_cmcs
from benchmark_tools.stats import save_stats, compute_stats

class AGSWrapper(Solver):


def main():

    params = ags_solver.Parameters()
    solver = ags_solver.Solver()
    solver.SetParameters(params)

    results = solve_class(grish_class(), AGSWrapper(solver))
    cmc_curve, avg_calculations = compute_stats(results)
    plot_cmc(cmc_curve)
    save_stats(cmc_curve, avg_calculations)

if __name__ == '__main__':
    main()
