export OMP_NUM_THREADS=28
./tpf_driver --seed 0 --type bootstrap --npart 40000 --output-file nkmp_p0_bootstrap.json
./tpf_driver --seed 0 --type resample --npart 40000 --output-file nkmp_p0_resample.json
./tpf_driver --seed 0 --npart 40000 --output-file nkmp_p0_tempered_40kr2.json
./tpf_driver --seed 0 --npart 40000 --rstar 3.0 --output-file nkmp_p0_tempered_40kr3.json
./tpf_driver --seed 0 --npart 4000 --output-file nkmp_p0_tempered_4kr2.json
./tpf_driver --seed 0 --npart 4000 --rstar 3.0 --output-file nkmp_p0_tempered_4kr3.json
./tpf_driver --seed 0 --type opt --npart 400 --output-file nkmp_p0_opt.json
python src/tpf_tables_and_figures.py nkmp_p0_bootstrap.json nkmp_p0_resample.json nkmp_p0_tempered_40kr2.json nkmp_p0_tempered_40kr3.json nkmp_p0_tempered_4kr2.json nkmp_p0_tempered_4kr3.json nkmp_p0_opt.json
#./tpf_driver --pmsv 1 --type bootstrap --npart 40000 --output-file nkmp_p0_bootstrap.json
#./tpf_driver --pmsv 1 --type resample --npart 40000 --output-file nkmp_p0_resample.json
#./tpf_driver --pmsv 1 --npart 40000 --output-file nkmp_p0_tempered_40kr2.json
#./tpf_driver --pmsv 1 --npart 40000 --rstar 3.0 --output-file nkmp_p0_tempered_40kr3.json
#./tpf_driver --pmsv 1 --npart 4000 --output-file nkmp_p0_tempered_4kr2.json
#./tpf_driver --pmsv 1 --npart 4000 --rstar 3.0 --output-file nkmp_p0_tempered_4kr3.json
#./tpf_driver --pmsv 1 --type opt --npart 400 --output-file nkmp_p0_opt.json
# python src/tpf_tables_and_figures.py nkmp_p0_bootstrap.json nkmp_p0_resample.json nkmp_p0_tempered_40kr2.json nkmp_p0_tempered_40kr3.json nkmp_p0_tempered_4kr2.json nkmp_p0_tempered_4kr3.json nkmp_p0_opt.json
