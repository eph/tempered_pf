#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
awk '/^processor/{n+=1}END{print n}' /proc/cpuinfo
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
#./tpf_driver --seed 0 --model sw --nsim 100 --type bootstrap --npart 40000 --output-file sw_p1_bootstrap.json
./tpf_driver --pmsv 1 --seed 0 --model sw --nsim 200 --npart 40000 --output-file sw_p1_tempered_40kr2.json
./tpf_driver --pmsv 1 --seed 0 --model sw --nsim 200 --npart 40000 --rstar 3.0 --output-file sw_p1_tempered_40kr3.json
./tpf_driver --pmsv 1 --seed 0 --model sw --nsim 200 --npart 6500 --output-file sw_p1_tempered_4kr2.json
./tpf_driver --pmsv 1 --seed 0 --model sw --nsim 200 --npart 7500 --rstar 3.0 --output-file sw_p1_tempered_4kr3.json
python src/tpf_tables_and_figures.py sw_p1_bootstrap.json sw_p1_tempered_4kr2.json sw_p1_tempered_4kr3.json sw_p1_tempered_40kr2.json sw_p1_tempered_40kr3.json > output/table5_part2.tex
