#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
awk '/^processor/{n+=1}END{print n}' /proc/cpuinfo
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
./tpf_driver --pmsv 1 --sample great_recession --seed 0 --nsim 200 --type bootstrap --npart 40000 --output-file nkmp_gr_p1_bootstrap.json
./tpf_driver --pmsv 1 --sample great_recession --seed 0 --nsim 200 --npart 40000 --output-file nkmp_gr_p1_tempered_40kr2.json
./tpf_driver --pmsv 1 --sample great_recession --seed 0 --nsim 200 --npart 40000 --rstar 3.0 --output-file nkmp_gr_p1_tempered_40kr3.json
./tpf_driver --pmsv 1 --sample great_recession --seed 0 --nsim 200 --npart 5500 --output-file nkmp_gr_p1_tempered_4kr2.json
./tpf_driver --pmsv 1 --sample great_recession --seed 0 --nsim 200 --npart 7000 --rstar 3.0 --output-file nkmp_gr_p1_tempered_4kr3.json
./tpf_driver --pmsv 1 --sample great_recession --seed 0 --nsim 200 --type resample --nintmh 10 --npart 40000 --output-file nkmp_gr_p1_resample.json
./tpf_driver --pmsv 1 --sample great_recession --seed 0 --nsim 200 --type opt --npart 400 --output-file nkmp_gr_p1_opt.json
python src/tpf_tables_and_figures.py nkmp_gr_p1_bootstrap.json nkmp_gr_p1_tempered_4kr2.json nkmp_gr_p1_tempered_4kr3.json nkmp_gr_p1_tempered_40kr2.json nkmp_gr_p1_tempered_40kr3.json  nkmp_gr_p1_resample.json nkmp_gr_p1_opt.json > output/table2_part2.tex
