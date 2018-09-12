import argparse
import json
import csv
import os
from tabulate import tabulate

def main(args):
    all_stats = []
    for subdir, dirs, files in sorted(os.walk(args.root_folder)):
        print('Prosessing directory ' + subdir)
        directory_stats = []
        for file in files:
            _, ext = os.path.splitext(file)
            if ext == '.json':
                full_path = os.path.join(subdir, file)
                stats = {}
                with open(full_path, 'r') as f:
                    stats = json.load(f)
                    directory_stats.append(stats)
        if len(directory_stats):
            all_stats.append((directory_stats, os.path.basename(subdir)))

    all_tables = {}
    columns = ['Метод', 'Среднее число испытаний', 'Решено задач']
    with open(os.path.join(args.root_folder, 'report.csv'), 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['class', 'method', 'solved', 'avg trials'])
        for class_results, class_name in all_stats:
            table_rows = []
            for method_results in class_results:
                table_rows.append((method_results['capture'], method_results['calc_counters'][0], method_results['num_solved']))
                writer.writerow([class_name, method_results['capture'],
                    method_results['num_solved'], round(method_results['calc_counters'][0])])
            writer.writerow([])
            table = tabulate(sorted(table_rows), headers=columns, tablefmt='latex', floatfmt='.2f', numalign='center')
            all_tables[class_name] = table

    all_lines = []
    with open(args.report_tex, 'r') as report_tex:
        all_lines = report_tex.readlines()
        for i, line in enumerate(all_lines):
            if '%table' in line:
                class_name = line.strip().split(' ')[-1]
                all_lines[i] = line.replace(line, line + all_tables[class_name])
        with open(os.path.join(os.path.split(args.report_tex)[0], 'report.tex'), 'w') as generated_report:
            for line in all_lines:
                generated_report.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('root_folder', type=str)
    parser.add_argument('--report_tex', type=str)
    main(parser.parse_args())
