# Dune Test	
This repository holds a basic DUNE module that solves a very basic Reaction Diffusion Problem of the two ODE's:

![ODE](images/ODE.png)

With parameters: 
	
* k1  = 1
* v1 = 0.06
* v2 = 0.12
* V = 0.5
* Km = 0.1

## Setup

Compilation of the project is possible on CYGWIN, OSX and Linux. As dependencies only `cmake` and `curl` are needed. From there simply run first: 

	./cloneAll.sh

which will checkout all Dune packages from their respective Dune repositories here on github. Next it might be possible to apply some patches: 

	cd dune-common && git apply ../patches/dune_common_ignore_fortran.patch && cd .. 
	cd dune-grid && git apply ../patches/dune_grid_cygwin.patch && cd ..
	cd dune-pdelab && git apply ../patches/dune_pdelab_new_cmake.patch && cd ..

After that all should be ready and 

	./compileDune.sh

should compile DUNE as well as three grid libraries (ALUgrid, Alberta, SuperLU) into `./install`. Next: 

	./compileTest.sh 

ought to compile the Dune-test module, and you should be able to run it by invoking 

	cd build-test/src
	./dune_test

## Acknowledgments
The model is based of a C++ implementation by Sven Sahle. The original Model has been developed by Pavel Hron.