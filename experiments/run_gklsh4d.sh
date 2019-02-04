declare -a algos=("ags")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../experiments/solve_different_methods.py --verbose --dist_stop --max_iters 250000 --problems_dim 4 --problems_class gklsh --algo $algo --stats_fname ../experiments/gklsh4d/$algo.json
done
#python2 ../experiments/solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture AGSl --verbose --max_iters 250000 --problems_dim 4 --problems_class gklsh --r 4.9 --lm 5 --eps 0.01 --stats_fname ../experiments/gklsh4d/agsl_examin.json
python2 ../experiments/plot_cmcs.py ../experiments/gklsh4d/ --show
