import sys
import re
import argparse
from subprocess import Popen, PIPE
from stats import save_stats, compute_stats

def create_parser ():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bin_path', type=str, default='C:\Projects\globalizer\_bin\examin_d.exe ')    
    parser.add_argument('--max_iters', type=int, default=10000, help='limit of iterations for the method')
    parser.add_argument('--problems_lib', type=str, default='gkls.dll')
    parser.add_argument('--problems_dim', type=int, default=2)
    parser.add_argument('--eps', type=float, default=0.01)
    parser.add_argument('--stopCond', type=str, default='OptimumVicinity')
    parser.add_argument('--stats_fname', type=str, default='test.json')
    parser.add_argument('--r', type=float, default='2.3')
    parser.add_argument('--lm', type=float, default='0')
    
    return parser

def create_start_command (namespace, function_number):
    string_start_examin = str(namespace.bin_path)
    string_start_examin += ' -MaxNP ' + str(namespace.max_iters)
    string_start_examin += ' -lib ' + str(namespace.problems_lib)
    string_start_examin += ' -function_number ' + str(function_number)
    string_start_examin += ' -N ' + str(namespace.problems_dim)
    string_start_examin += ' -Eps ' + str(namespace.eps)
    string_start_examin += ' -stopCond ' + str(namespace.stopCond)
    string_start_examin += ' -r ' + str(namespace.r)
    string_start_examin += ' -lm '+ str(namespace.lm)
    
    return string_start_examin

def start_examin (args, function_number):
    string_start_examin = create_start_command(args, function_number)
    proc = Popen(string_start_examin, shell=True, stdout=PIPE, stderr=PIPE)
    proc.wait()
    output_examin = proc.communicate()

    return output_examin[0]

def start_serial(args):
    information_about_experiments = []
    problem_status = []
    for i in range(1,101):
        try:
            output_examin = start_examin(args, i)
        
            number_of_trials = re.search(r'NumberOfTrials = (\d+)', str(output_examin)).group(1)
            information_about_experiments.append([int(number_of_trials)])
            result = re.search(r'FOUND!', str(output_examin))
            if (result==None):
                problem_status.append(0)
            else:
                problem_status.append(1)
            print(i)
        except BaseException:
            print(output_examin)
            sys.exit()
    stats = compute_stats(information_about_experiments, problem_status)

    return stats
 
def main(args):
    stats = start_serial(args)
    save_stats(stats, args.stats_fname)

if __name__ == '__main__':
    parser = create_parser()

    main(parser.parse_args())
