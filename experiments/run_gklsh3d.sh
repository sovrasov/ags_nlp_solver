declare -a algos=("ags")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../experiments/solve_different_methods.py --verbose --dist_stop --max_iters 25000 --problems_dim 3 --problems_class gklsh --algo $algo --stats_fname ../experiments/gklsh3d/$algo.json
done
#python2 ../experiments/solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture AGSl --verbose --max_iters 25000 --problems_dim 3 --problems_class gklsh --r 4.4 --lm 5 --eps 0.01 --stats_fname ../experiments/gklsh3d/agsl_examin.json
python2 ../experiments/plot_cmcs.py ../experiments/gklsh3d/ --show
