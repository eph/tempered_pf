#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
awk '/^processor/{n+=1}END{print n}' /proc/cpuinfo
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
#./tpf_driver --seed 0 --pmsv 1 --model sw --nsim 100 --type bootstrap --npart 40000 --output-file sw_p1_bootstrap.json
./tpf_driver --seed 0 --pmsv 1 --model sw --nsim 100 --npart 40000 --output-file sw_p1_tempered_40kr2_nmh10.json
./tpf_driver --seed 0 --pmsv 1 --model sw --nsim 100 --npart 40000 --rstar 3.0 --output-file sw_p1_tempered_40kr3_nmh10.json
./tpf_driver --seed 0 --pmsv 1 --model sw --nsim 100 --npart 2000 --output-file sw_p1_tempered_4kr2_nmh10.json
./tpf_driver --seed 0 --pmsv 1 --model sw --nsim 100 --npart 2000 --rstar 3.0 --output-file sw_p1_tempered_4kr3_nmh10.json
./tpf_driver --seed 0 --pmsv 1 --model sw --type resample --nintmh 10 --npart 40000 --output-file sw_p1_resample.json
./tpf_driver --seed 0 --pmsv 1 --model sw --type opt --npart 400 --output-file sw_p1_opt.json
python src/tpf_tables_and_figures.py sw_p1_bootstrap.json sw_p1_tempered_4kr2_nmh10.json sw_p1_tempered_4kr3_nmh10.json sw_p1_tempered_40kr2_nmh10.json sw_p1_tempered_40kr3_nmh10.json  sw_p1_resample.json sw_p1_opt.json > table6_part2.tex
