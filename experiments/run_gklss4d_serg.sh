declare -a algos=("ags" "direct" "directl" "mlsl" "crs" "scd" "stogo", "sda")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../experiments/solve_different_methods.py --serg_eps --verbose --dist_stop --max_iters 150000 --problems_dim 4 --problems_class gklss --algo $algo --stats_fname ../experiments/gklss4d_serg/$algo.json
done
python2 ../experiments/solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture AGSl --verbose --max_iters 150000 --problems_dim 4 --problems_class gklss --r 4.7 --lm 5 --serg_eps --stats_fname ../experiments/gklss4d_serg/agsl_examin.json
python2 ../experiments/plot_cmcs.py ../experiments/gklss4d_serg/ --show
