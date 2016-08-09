make clean
make tpf_driver
cp tpf_driver $PREFIX/bin
mkdir -p $PREFIX/include/tempered_pf
mkdir -p $PREFIX/include/tempered_pf/nkmp
cp src/nkmp/us.txt $PREFIX/include/tempered_pf/nkmp
cp p0.txt $PREFIX/include/tempered_pf
