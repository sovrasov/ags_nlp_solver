declare -a algos=("ags") # "direct" "directl" "mlsl" "crs" "scd" "stogo", "simple", "sda")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../experiments/solve_different_methods.py --verbose --dist_stop --max_iters 9000 --problems_dim 2 --problems_class gklsh --algo $algo --stats_fname ../experiments/gklsh2d/ags_r.json
done
#python2 ../experiments/solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture AGSl --verbose --max_iters 9000 --problems_dim 2 --problems_class gklsh --r 6.5 --lm 5 --eps 0.01 --stats_fname ../experiments/gklsh2d/agsl_examin.json
#python2 ../experiments/plot_cmcs.py ../experiments/gklsh2d/ --show
