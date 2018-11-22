#!/bin/bash

declare -a n_nodes=("1" "2" "4" "8" "12")
r_values=( ["1"]="4.7" ["2"]="4.4" ["4"]="4.1" ["8"]="3.8" ["12"]="3.5")
declare -a thr_nums=("1" "16")

iters=150000
class="gklss"
folder="mpi/gklss4d"
partition="-p gpu --reservation=global_opt"
complexity="20000"

for n_num in "${n_nodes[@]}"
do
preffix="srun -N $n_num $partition -t 600"
for thr_num in "${thr_nums[@]}"
  do
    python2 solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture 'Parallel AGS' --verbose --max_iters $iters --eps 0.01 --problems_class $class --problems_dim 4 --lm 0 --r "${r_values[$n_num]}" --stats_fname $folder/examin_"$n_num"_"$thr_num".json --parallel_evolvent --ne $n_num --p $thr_num --preffix ''"$preffix"'' --complexity $complexity
  done
done
