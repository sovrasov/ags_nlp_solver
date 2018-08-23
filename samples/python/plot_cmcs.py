import argparse
import glob
import os
from benchmark_tools.plot import plot_cmcs
from benchmark_tools.stats import load_stats, compute_stats

def main(args):
    source_files = glob.glob(os.path.join(args.stats_folder, '*.json'))
    all_stats = [load_stats(file) for file in source_files]
    captures = [item['capture'] for item in all_stats]
    iters = [item['calc_counters'][-1] for item in all_stats]
    cmcs = [[item['cmc_iters'], item['cmc_vals']] for item in all_stats]
    print(zip(captures, iters))

    plot_cmcs(cmcs, captures=captures, show=args.show,
              filename=args.output if len(args.output) > 0 else os.path.join(args.stats_folder, 'cmc.pdf'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('stats_folder', type=str)
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--output', type=str, default='')

    main(parser.parse_args())
