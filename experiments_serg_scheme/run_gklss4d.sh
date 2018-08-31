declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../samples/python/solve_different_methods.py --serg_eps --verbose --dist_stop --max_iters 150000 --problems_dim 4 --problems_class gklss --algo $algo --stats_fname ../experiments_serg_scheme/gklss4d/$algo.json
done
python2 ../samples/python/plot_cmcs.py ../experiments_serg_scheme/gklss4d/ --show
