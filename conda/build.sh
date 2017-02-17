make clean
rm -f tpf_driver
make tpf_driver
cp tpf_driver $PREFIX/bin
#make test
#make test_driver
#cp test_driver $PREFIX/bin/tpf_test
cat tpf_driver_tmp.f90
#./test_driver
mkdir -p $PREFIX/include/tempered_pf
mkdir -p $PREFIX/include/tempered_pf/nkmp
cp src/nkmp/us.txt $PREFIX/include/tempered_pf/nkmp
cp p0.txt $PREFIX/include/tempered_pf/nkmp
