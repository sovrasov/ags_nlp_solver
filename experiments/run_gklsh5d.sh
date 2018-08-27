declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../samples/python/solve_different_methods.py --verbose --dist_stop --max_iters 600000 --problems_dim 5 --problems_class gklsh --algo $algo --stats_fname ../experiments/gklsh5d/$algo.json
done
