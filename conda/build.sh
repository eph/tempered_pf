make clean
make tpf_driver_nkmp
cp tpf_driver_nkmp $PREFIX/bin
mkdir -p $PREFIX/include/tempered_pf
cp src/us.txt $PREFIX/include/tempered_pf
cp p0.txt $PREFIX/include/tempered_pf
