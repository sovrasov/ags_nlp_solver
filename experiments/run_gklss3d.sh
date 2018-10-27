declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo", "simple", "sda")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../experiments/solve_different_methods.py --verbose --dist_stop --max_iters 15000 --problems_dim 3 --problems_class gklss --algo $algo --stats_fname ../experiments/gklss3d/$algo.json
done
python2 ../experiments/solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture AGSl --verbose --max_iters 15000 --problems_dim 3 --problems_class gklss --r 3.7 --lm 5 --eps 0.01 --stats_fname ../experiments/gklss3d/agsl_examin.json
python2 ../experiments/plot_cmcs.py ../experiments/gklss3d/ --show
