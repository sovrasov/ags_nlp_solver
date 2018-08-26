declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../samples/python/solve_different_methods.py --verbose --dist_stop --max_iters 5000 --problems_dim 2 --problems_class grish --algo $algo #--stats_fname ../experiments/grish/$algo.json
done
