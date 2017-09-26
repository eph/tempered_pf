#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
awk '/^processor/{n+=1}END{print n}' /proc/cpuinfo
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
./tpf_driver --seed 0  --model sw --nsim 200 --type bootstrap --npart 40000 --output-file sw_p0_bootstrap.json
./tpf_driver --seed 0  --model sw --nsim 200 --npart 2000 --nintmh 10 --output-file sw_p0_tempered_4kr2_nmh10.json
./tpf_driver --seed 0  --model sw --nsim 200 --npart 2800 --nintmh 10 --rstar 3.0 --output-file sw_p0_tempered_4kr3_nmh10.json
./tpf_driver --seed 0  --model sw --nsim 200 --npart 40000 --nintmh 10 --output-file sw_p0_tempered_40kr2_nmh10.json
./tpf_driver --seed 0  --model sw --nsim 200 --npart 40000 --nintmh 10 --rstar 3.0 --output-file sw_p0_tempered_40kr3_nmh10.json
./tpf_driver --seed 0  --model sw --nsim 200 --type resample --nintmh 10 --npart 40000 --output-file sw_p0_resample.json
./tpf_driver --seed 0  --model sw --nsim 200 --type opt --npart 400 --output-file sw_p0_opt.json
python src/tpf_tables_and_figures.py sw_p0_bootstrap.json sw_p0_tempered_4kr2_nmh10.json sw_p0_tempered_4kr3_nmh10.json sw_p0_tempered_40kr2_nmh10.json sw_p0_tempered_40kr3_nmh10.json  sw_p0_resample.json sw_p0_opt.json > output/table6_part1.tex
