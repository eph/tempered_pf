SRC=src
TEST=tests
VPATH=.:$(SRC):$(TEST)

FC=mpif90

LOBJS=model_nkmp.o
ifdef CONDA_BUILD
LIB=$(PREFIX)/lib
INC=$(PREFIX)/include
else
LIB=$(HOME)/anaconda3/lib
INC=$(HOME)/anaconda3/include
endif

FPP=fypp
FRUIT=-I$(INC)/fruit -L$(LIB) -lfruit
FLAP=-I$(INC)/flap -L$(LIB) -lflap
FORTRESS=-I/home/eherbst/Dropbox/code/fortress -L/home/eherbst/Dropbox/code/fortress -lfortress

%.o : %.f90
	$(FPP) -DGFORTRAN $< $(notdir $(basename $<))_tmp.f90
	$(FC) $(FRUIT) $(FLAP) $(FORTRESS) -fPIC -c $(notdir $(basename $<)_tmp.f90) -o $(notdir $(basename $<)).o
	rm $(notdir $(basename $<))_tmp.f90

test_driver: test_driver.f90 $(LOBJS) test_nkmp.o
	$(FC) $^  $(FRUIT) $(FORTRESS) $(FLAP) -llapack  -o $@ 

test:
	python conda/run_test.py

clean:
	rm -f *.o *.mod
