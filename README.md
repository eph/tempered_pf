Companion Code for _Tempered Particle Filtering_
================================================
by Ed Herbst and Frank Schorfheide

Click [here](#) to read the latest version of the draft.

Installation 
------------ 

If you're using a linux distribution, the easiest way to install this
software is by using [Conda](http://www.continuum.io/downloads), a
packaging tool which helps disseminate scientific software.

At the command prompt, type: 

```
conda --add channels conda-forge eherbst 
conda install tempered_pf
```

Conda will install executables `tpf_everything` `tpf_driver`, `tpf_figures_and_tables`.

### Installation by hand

This project is written principally in Fortran, and so requires a
fortran compiler.  It uses the `fortress` library (available
[here](http://github.com/eph/fortress).)  Installation goes like:

1. Install `fortress` by hand or via Conda.
2. Clone / download this repository.
3. Edit the `makefile` to link to libfortress.so (and its dependencies) correctly. 
4. At the prompt:
```
make tpf
```

This will result in the 4 executables being placed in the `bin/` directory.


Usage 
----- 
The main program is tpf driver, which runs all of the calculations
reported in the paper.  

```
eherbst@thnkpd:~$ tpf_driver --help
usage: Tempered Particle Filtering Example  [--model value] [--sample value] [--npart value] [--pmsv value] [--nintmh value] [--rstar value] [--nsim value] [--output-file value] [--help] [--version]

Optional switches:
   --model value, -m value, value in: `nkmp,sw`
    default value nkmp
    Model
   --sample value, -s value, value in: `great_moderation,great_recession`
    default value great_moderation
    Sample
   --npart value, -n value
    default value 4000
    Number of particles
   --pmsv value, -p0 value
    default value p0.txt
    Parameter File
   --nintmh value, -mh value
    default value 1
    Number of intermediate MH steps
   --rstar value, -r value
    default value 2.0
    Ineffeciency Ration (for TPF)
   --nsim value
    default value 100
    Ineffeciency Ration (for TPF)
   --output-file value, -o value
    default value output.json
    Output File
   --help, -h
    Print this help message
   --version, -v
    Print version
```


