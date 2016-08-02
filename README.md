Companion Code for _Tempered Particle Filtering_
================================================
by Ed Herbst and Frank Schorfheide

Click [here](#) to read the latest version of the draft.

Installation
------------
If you're using a linux distribution, the easiest way to install software is using [Conda](http://www.continuum.io/downloads), a packaging tool which helps disseminate scientific software. 

At the command prompt, type: 

```
conda install -c eherbst tempered_pf
```

Conda will install executables `tpf_everything` `tpf_driver_nkmp`,
`tpf_driver_sw`, `tpf_figures_and_tables`.

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
usage instructions here


