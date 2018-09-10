declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo", "sda")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../samples/python/solve_different_methods.py --serg_eps --verbose --dist_stop --max_iters 250000 --problems_dim 4 --problems_class gklsh --algo $algo --stats_fname ../experiments/gklsh4d_serg/$algo.json
done
python2 ../samples/python/plot_cmcs.py ../experiments/gklsh4d_serg/ --show
