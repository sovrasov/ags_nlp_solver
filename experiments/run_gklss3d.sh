declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo", "simple")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../samples/python/solve_different_methods.py --verbose --dist_stop --max_iters 15000 --problems_dim 3 --problems_class gklss --algo $algo --stats_fname ../experiments/gklss3d/$algo.json
done
python2 ../samples/python/plot_cmcs.py ../experiments/gklss3d/ --show
