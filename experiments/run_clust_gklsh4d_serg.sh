#!/bin/bash

declare -a n_nodes=("1" "2" "4" "8" "12")
r_values=( ["1"]="4.9" ["2"]="4.6" ["4"]="4.4" ["8"]="4.0" ["12"]="3.7")
declare -a thr_nums=("1" "16")

iters=250000
class="gklsh"
folder="mpi/gklsh4d_serg"
partition="-p gpu --reservation=global_opt"
complexity="20000"

for n_num in "${n_nodes[@]}"
do
preffix="srun -N $n_num $partition -t 600"
for thr_num in "${thr_nums[@]}"
  do
    python2 solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture 'Parallel AGS' --verbose --max_iters $iters --serg_eps --problems_class $class --problems_dim 4 --lm 0 --r "${r_values[$n_num]}" --stats_fname $folder/examin_"$n_num"_"$thr_num".json --parallel_evolvent --ne $n_num --p $thr_num --preffix ''"$preffix"'' --complexity $complexity
  done
done
