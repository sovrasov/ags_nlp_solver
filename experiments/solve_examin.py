import sys, os, math
import re
import argparse
from subprocess import Popen, PIPE
from benchmark_tools.stats import save_stats, compute_stats
import platform

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bin_path', type=str, default='')

    parser.add_argument('--problems_class', type=str, choices=['grish','gklss','gklsh'], default='grish')
    parser.add_argument('--problems_dim', type=int, default=2)
    parser.add_argument('--stop_cond', type=str, default='OptimumVicinity')
    parser.add_argument('--serg_eps', action='store_true')

    parser.add_argument('--max_iters', type=int, default=10000, help='limit of iterations for the method')
    parser.add_argument('--eps', type=float, default=0.01)
    parser.add_argument('--r', type=float, default='2.3')
    parser.add_argument('--lm', type=int, default='0')
    parser.add_argument('--m', type=int, default='10')
    parser.add_argument('--p', type=int, default='1', help='Number of OMP threads')
    parser.add_argument('--ne', type=int, default='1', help='Number of evolvents')
    parser.add_argument('--stats_fname', type=str, default='test.json')
    parser.add_argument('--verbose', action='store_true', help='Print additional info to console')
    parser.add_argument('--parallel_evolvent', action='store_true', help='Enable parallelism by evolvents')
    parser.add_argument('--algo_capture', type=str, default='AGS-Examin')
    parser.add_argument('--preffix', type=str, default='')

    return parser

def create_parameters_dict(cl_args):
    parameters = {}
    parameters['MaxNP'] = cl_args.max_iters
    parameters['lib'] = os.path.join(os.path.abspath(cl_args.bin_path), \
        get_platform_lib_name(get_lib_name_by_functions_class(cl_args.problems_class)))
    parameters['function_number'] = 1
    parameters['N'] = cl_args.problems_dim
    parameters['stopCond'] = cl_args.stop_cond
    parameters['r'] = cl_args.r
    parameters['lm'] = cl_args.lm
    parameters['m'] = cl_args.m
    parameters['spm'] = cl_args.max_iters // 10
    parameters['nt'] = cl_args.p
    parameters['np'] = cl_args.p
    parameters['ml'] = cl_args.ne
    if cl_args.ne > 1:
        parameters['tm'] = 'ParallelMultievolventsMethod'
        parameters['mt'] = 'mpRotated'
        if cl_args.parallel_evolvent:
            parameters['mpl'] = cl_args.ne

    if cl_args.serg_eps:
        serg_eps = {2: 0.01, 3: 0.01, 4: math.pow(1e-6, 1./4), 5: math.pow(1e-7, 1./5)}
        parameters['Eps'] = serg_eps[cl_args.problems_dim]
    else:
        parameters['Eps'] = cl_args.eps

    if cl_args.problems_class == 'gklss':
        parameters['PC'] = 'Simple'
    elif cl_args.problems_class == 'gklsh':
        parameters['PC'] = 'Hard'

    return parameters

def create_start_command(bin_path, parameters_dict):
    string_start_examin = os.path.join(os.path.abspath(bin_path), get_platform_executable_name('examin'))

    for param in parameters_dict:
        string_start_examin += ' -' + param + ' ' + str(parameters_dict[param])

    return string_start_examin

def start_examin(preffix, bin_path, parameters_dict):
    string_start_examin = preffix + ' ' + create_start_command(bin_path, parameters_dict)
    proc = Popen(string_start_examin, shell=True, stdout=PIPE, stderr=PIPE)
    proc.wait()
    output_examin = proc.communicate()
    if len(output_examin[0]) == 0:
        raise Exception('Examin launch failed: ' + string_start_examin)

    return output_examin[0]

def get_lib_name_by_functions_class(class_name):
    if 'grish' in class_name:
        return 'grishagin'
    elif 'gkls' in class_name:
        return 'gkls'
    else:
        raise Exception('Invalid functions class name')

def get_platform_lib_name(lib):
    if 'Linux' in platform.system():
        return 'lib' + lib + '.so'
    elif 'Windows' in platform.system():
        return lib + '.dll'
    return None

def get_platform_executable_name(name):
    if 'Linux' in platform.system():
        return name
    elif 'Windows' in platform.system():
        return name + '.exe'
    return None

def start_serial(args):
    information_about_experiments = []
    exec_times = []
    problem_status = []
    parameters = create_parameters_dict(args)
    print('Launch parameters: {}'.format(parameters))

    for i in range(1,101):
        output_examin = ''
        try:
            parameters['function_number'] = i
            output_examin = str(start_examin(args.preffix, args.bin_path, parameters))

            number_of_trials = re.search(r'NumberOfTrials = (\d+)', output_examin).group(1)
            information_about_experiments.append([int(number_of_trials)])
            time = float(re.search(r'Solve time = (\d+.\d+)', output_examin).group(1))
            exec_times.append(time)
            result = re.search(r'FOUND!', str(output_examin))
            if result == None:
                problem_status.append(False)
            else:
                problem_status.append(True)

            if args.verbose:
                print('Problem # {}: '.format(i) + ('solved' if problem_status[-1] else 'not solved') +
                    ' after {} trials'.format(number_of_trials))
        except BaseException as e:
            print(e)
            print(output_examin)
            sys.exit()
    stats = compute_stats(information_about_experiments, problem_status, exec_times)
    print('Problems solved: {}'.format(stats['num_solved']))
    for i, avg in enumerate(stats['avg_calcs'][:-1]):
        print('Average number of calculations of constraint #{}: {}'.format(i, avg))
    print('Average number of calculations of objective: {}'.format(stats['avg_calcs'][-1]))
    print('Average solver run time {}'.format(stats['avg_time']))

    return stats

def main(args):
    stats = start_serial(args)
    save_stats(stats, args.stats_fname, args.algo_capture)

if __name__ == '__main__':
    parser = create_parser()

    main(parser.parse_args())
