make clean
make tpf_driver
cp tpf_driver $PREFIX/bin
mkdir -p $PREFIX/include/tempered_pf
cp src/us.txt $PREFIX/include/tempered_pf
cp p0.txt $PREFIX/include/tempered_pf
