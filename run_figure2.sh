#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
export OMP_NUM_THREADS=20
./tpf_driver --seed 23 --nsim 100 --type bootstrap --npart 40000 --output-file nkmp_p0_bootstrap.json
./tpf_driver --seed 23 --nsim 10 --npart 40000 --output-file nkmp_p0_tempered_40kr2.json
# ./tpf_driver --seed 23 --nsim 100 --save-states --npart 40000 --rstar 3.0 --output-file nkmp_p0_tempered_40kr3.json
# ./tpf_driver --seed 23 --nsim 100 --save-states --npart 4000 --output-file nkmp_p0_tempered_4kr2.json
# ./tpf_driver --seed 23 --nsim 100 --save-states --npart 4000 --rstar 3.0 --output-file nkmp_p0_tempered_4kr3.json
# ./tpf_driver --seed 23 --nsim 100 --type resample --nintmh 10 --npart 40000 --output-file nkmp_p0_resample.json
# ./tpf_driver --seed 23 --type opt --npart 400 --output-file nkmp_p0_opt.json
#python src/tpf_tables_and_figures.py nkmp_p0_bootstrap.json nkmp_p0_tempered_4kr2.json nkmp_p0_tempered_4kr3.json nkmp_p0_tempered_40kr2.json nkmp_p0_tempered_40kr3.json  nkmp_p0_resample.json nkmp_p0_opt.json > table2_part1.tex
#python src/tpf_std_filtered_state.py nkmp_p0_bootstrap.json nkmp_p0_tempered_40kr2.json nkmp_p0_tempered_40kr3.json nkmp_p0_tempered_4kr2.json nkmp_p0_tempered_4kr3.json 
# ./tpf_driver --pmsv 1 --seed 0 --nsim 100 --save-states --type bootstrap --npart 40000 --output-file nkmp_p0_bootstrap.json
# ./tpf_driver --pmsv 1 --seed 0 --nsim 100 --save-states --npart 40000 --output-file nkmp_p0_tempered_40kr2.json
# ./tpf_driver --pmsv 1 --seed 0 --nsim 100 --save-states --npart 40000 --rstar 3.0 --output-file nkmp_p0_tempered_40kr3.json
# ./tpf_driver --pmsv 1 --seed 0 --nsim 100 --save-states --npart 4000 --output-file nkmp_p0_tempered_4kr2.json
# ./tpf_driver --pmsv 1 --seed 0 --nsim 100 --save-states --npart 4000 --rstar 3.0 --output-file nkmp_p0_tempered_4kr3.json
# ./tpf_driver --pmsv 1 --seed 0 --type resample --npart 40000 --output-file nkmp_p0_resample.json
# ./tpf_driver --pmsv 1 --seed 0 --type opt --npart 400 --output-file nkmp_p0_opt.json
# python src/tpf_tables_and_figures.py nkmp_p0_bootstrap.json nkmp_p0_tempered_4kr2.json nkmp_p0_tempered_4kr3.json nkmp_p0_tempered_40kr2.json nkmp_p0_tempered_40kr3.json  nkmp_p0_resample.json nkmp_p0_opt.json > table2_part2.tex


# Great Recession Sample
# ./tpf_driver --sample great_recession --seed 0 --nsim 100 --type bootstrap --npart 40000 --output-file nkmp_p0_bootstrap.json
# ./tpf_driver --sample great_recession --seed 0 --nsim 100 --npart 40000 --output-file nkmp_p0_tempered_40kr2.json
# ./tpf_driver --sample great_recession --seed 0 --nsim 100 --npart 40000 --rstar 3.0 --output-file nkmp_p0_tempered_40kr3.json
# ./tpf_driver --sample great_recession --seed 0 --nsim 100 --npart 4000 --output-file nkmp_p0_tempered_4kr2.json
# ./tpf_driver --sample great_recession --seed 0 --nsim 100 --npart 4000 --rstar 3.0 --output-file nkmp_p0_tempered_4kr3.json
# ./tpf_driver --sample great_recession --seed 0 --pmsv 0 --type resample --npart 40000 --output-file nkmp_p0_resample.json
# ./tpf_driver --sample great_recession --seed 0 --pmsv 0 --type opt --npart 400 --output-file nkmp_p0_opt.json
# python src/tpf_tables_and_figures.py nkmp_p0_bootstrap.json nkmp_p0_tempered_4kr2.json nkmp_p0_tempered_4kr3.json nkmp_p0_tempered_40kr2.json nkmp_p0_tempered_40kr3.json  nkmp_p0_resample.json nkmp_p0_opt.json > table3_part1.tex
# python src/tpf_evolution_figures.py nkmp_p0_tempered_40kr2.json --output evolution.pdf
# ./tpf_driver --sample great_recession --pmsv 1 --seed 0 --nsim 100 --type bootstrap --npart 40000 --output-file nkmp_p0_bootstrap.json
# ./tpf_driver --sample great_recession --pmsv 1 --seed 0 --nsim 100 --npart 40000 --output-file nkmp_p0_tempered_40kr2.json
# ./tpf_driver --sample great_recession --pmsv 1 --seed 0 --nsim 100 --npart 40000 --rstar 3.0 --output-file nkmp_p0_tempered_40kr3.json
# ./tpf_driver --sample great_recession --pmsv 1 --seed 0 --nsim 100 --npart 4000 --output-file nkmp_p0_tempered_4kr2.json
# ./tpf_driver --sample great_recession --pmsv 1 --seed 0 --nsim 100 --npart 4000 --rstar 3.0 --output-file nkmp_p0_tempered_4kr3.json
# ./tpf_driver --sample great_recession --pmsv 1 --seed 0 --type resample --npart 40000 --output-file nkmp_p0_resample.json
# ./tpf_driver --sample great_recession --pmsv 1 --seed 0 --type opt --npart 400 --output-file nkmp_p0_opt.json
# python src/tpf_tables_and_figures.py nkmp_p0_bootstrap.json nkmp_p0_tempered_4kr2.json nkmp_p0_tempered_4kr3.json nkmp_p0_tempered_40kr2.json nkmp_p0_tempered_40kr3.json  nkmp_p0_resample.json nkmp_p0_opt.json > table3_part2.tex
