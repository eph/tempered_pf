make clean
rm -f tpf_driver
make tpf_driver
cp tpf_driver $PREFIX/bin
mkdir -p $PREFIX/include/tempered_pf
mkdir -p $PREFIX/include/tempered_pf/nkmp
cp include/tempered_pf/nkmp/* $PREFIX/include/tempered_pf/nkmp

