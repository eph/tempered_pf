#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
awk '/^processor/{n+=1}END{print n}' /proc/cpuinfo
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
./tpf_driver  --seed 0 --nsim 200 --type bootstrap --npart 40000 --output-file nkmp_p0_bootstrap.json
./tpf_driver  --seed 0 --nsim 200 --npart 40000 --output-file nkmp_p0_tempered_40kr2.json
./tpf_driver  --seed 0 --nsim 200 --npart 40000 --rstar 3.0 --output-file nkmp_p0_tempered_40kr3.json
./tpf_driver  --seed 0 --nsim 200 --npart 7000 --output-file nkmp_p0_tempered_4kr2.json
./tpf_driver  --seed 0 --nsim 200 --npart 8500 --rstar 3.0 --output-file nkmp_p0_tempered_4kr3.json
./tpf_driver  --seed 0 --nsim 200 --type resample --nintmh 10 --npart 40000 --output-file nkmp_p0_resample.json
./tpf_driver  --seed 0 --nsim 200 --type opt --npart 400 --output-file nkmp_p0_opt.json
python src/tpf_tables_and_figures.py nkmp_p0_bootstrap.json nkmp_p0_tempered_4kr2.json nkmp_p0_tempered_4kr3.json nkmp_p0_tempered_40kr2.json nkmp_p0_tempered_40kr3.json  nkmp_p0_resample.json nkmp_p0_opt.json > output/table1_part1.tex
