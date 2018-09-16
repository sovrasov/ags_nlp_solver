declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo", "sda")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../experiments/solve_different_methods.py --verbose --dist_stop --max_iters 25000 --problems_dim 3 --problems_class gklsh --algo $algo --stats_fname ../experiments/gklsh3d/$algo.json
done
python2 ../experiments/plot_cmcs.py ../experiments/gklsh3d/ --show
