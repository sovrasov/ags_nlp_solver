declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../samples/python/solve_different_methods.py --verbose --dist_stop --max_iters 150000 --problems_dim 4 --problems_class gklss --algo $algo --stats_fname ../experiments/gklss4d/$algo.json
done
