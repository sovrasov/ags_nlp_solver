declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo", "simple")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../samples/python/solve_different_methods.py --verbose --dist_stop --max_iters 9000 --problems_dim 2 --problems_class gklsh --algo $algo --stats_fname ../experiments/gklsh2d/$algo.json
done
python2 ../samples/python/plot_cmcs.py ../experiments/gklsh2d/ --show
