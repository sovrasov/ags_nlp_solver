declare -a algos=("ags")

for algo in "${algos[@]}"
do
   echo "$algo"
   python2 ../experiments/solve_different_methods.py --serg_eps --verbose --dist_stop --max_iters 350000 --problems_dim 5 --problems_class gklss --algo $algo --stats_fname ../experiments/gklss5d_serg/$algo.json
done
#python2 ../experiments/solve_examin.py --bin_path ../../globalizer/_bin/ --algo_capture AGSl --verbose --max_iters 350000 --problems_dim 5 --problems_class gklss --r 4 --lm 5 --serg_eps --stats_fname ../experiments/gklss5d_serg/agsl_examin.json
python2 ../experiments/plot_cmcs.py ../experiments/gklss5d_serg/ --show
