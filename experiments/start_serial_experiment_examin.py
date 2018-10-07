import subprocess
import re
import argparse
from subprocess import Popen, PIPE
from stats import save_stats, compute_stats

def createParser ():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bin_path', type=str, default='C:\Projects\globalizer\_bin\examin_d.exe ')    
    parser.add_argument('--max_iters', type=int, default=10000, help='limit of iterations for the method')
    parser.add_argument('--problems_lib', type=str, default='gkls.dll')
    parser.add_argument('--function_number', type=int, default='1')
    parser.add_argument('--problems_dim', type=int, default=2)
    parser.add_argument('--eps', type=float, default=0.01)
    parser.add_argument('--stopCond', type=str, default='OptimumVicinity')
    parser.add_argument('--stats_fname', type=str, default='test')
    parser.add_argument('--r', type=float, default='2.3')
    parser.add_argument('--lm', type=float, default='0')
    
    return parser

def createStartCommand (namespace):
    string_start_examin = str(namespace.bin_path)
    string_start_examin += " -MaxNP " + str(namespace.max_iters)
    string_start_examin += " -lib " + str(namespace.problems_lib)
    string_start_examin += " -function_number " + str(namespace.function_number)
    string_start_examin += " -N " + str(namespace.problems_dim)
    string_start_examin += " -Eps " + str(namespace.eps)
    string_start_examin += " -stopCond " + str(namespace.stopCond)
    string_start_examin += " -r " + str(namespace.r)
    string_start_examin += " -lm " + str(namespace.lm)
    
    return string_start_examin

def startExamin (args):
    string_start_examin = createStartCommand(args)
    proc = Popen(string_start_examin, shell=True, stdout=PIPE, stderr=PIPE)
    proc.wait()
    outputExamin = proc.communicate()

    return outputExamin[0]

def startSerial(args):
    informationAboutExperiments = []
    problemStatus = []
    for i in range(1,101):
        args.function_number = i
        outputExamin = startExamin(args)
        
        NumberOfTrials = re.search(r'NumberOfTrials = (\d+)', str(outputExamin)).group(1)
        informationAboutExperiments.append([int(NumberOfTrials)])
        result = re.search(r'FOUND!', str(outputExamin))
        if (result==None):
            problemStatus.append(0)
        else:
            problemStatus.append(1)
        print(i)
    stats = compute_stats(informationAboutExperiments, problemStatus)

    return stats
 
def main(args):
    stats = startSerial(args)
    save_stats(stats, args.stats_fname)

if __name__ == '__main__':
    parser = createParser()

    main(parser.parse_args())
