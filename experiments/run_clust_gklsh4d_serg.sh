#!/bin/bash

declare -a n_nodes=("1" "2" "4" "8" "12")
declare -a thr_nums=("1" "16")

iters=250000
class="gklsh"
folder="mpi/gklsh4d_serg"
patrition="-p cpu"
complexity="0"

for n_num in "${n_nodes[@]}"
do
for thr_num in "${thr_nums[@]}"
  do
    python2 solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture 'Parallel AGS' --verbose --max_iters $iters --serg_eps --problems_class $class --problems_dim 4 --lm 0 --r 4.9 --stats_fname $folder/examin_"$n_num"_"$thr_num".json --parallel_evolvent --ne $n_num --p $thr_num --preffix 'srun -N $n_num $patrition' --complexity $complexity
  done
done
