declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo", "simple", "sda")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../experiments/solve_different_methods.py --verbose --dist_stop --max_iters 5000 --problems_dim 2 --problems_class grish --algo $algo --stats_fname ../experiments/grish/$algo.json
done
python2 ../experiments/plot_cmcs.py ../experiments/grish/ --show
