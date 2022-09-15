# Tic-tac: Two is company, three's a crowd
Research code for the simulation of three-nucleon scattering in preparation for open-source use.

## Summary
Tic-tac uses the [wave-packet continuum discretization](https://www.sciencedirect.com/science/article/abs/pii/S0003491615001773) (WPCD) method to solve the three-nucleon [Alt-Grassberger-Sandhas](https://www.sciencedirect.com/science/article/abs/pii/0550321367900168) (AGS) equations for elastic nucleon-deuteron scattering.
Currently, Tic-tac can produce on-shell neutron-deuteron scattering amplitudes (U-matrix elements).

## Table of contents
 - Summary
 - Dependencies
 - Compiling
 - Running the code
 - How to cite

## Dependencies
Tic-tac depends on the following libraries:
 - [BLAS](https://netlib.org/blas/)
 - [LAPACK](https://netlib.org/lapack/)
 - [OpenMP](https://www.openmp.org/)
 - [GSL](https://www.gnu.org/software/gsl/)
 - [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
 
Tic-tac relies on the user parsing nuclear-potential codes directly into the source code.
It comes equipped with several codes for different nuclear potential, notably among which are:
 - [Nijmegen-I](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.49.2950)
 - [N2LO<sub>opt](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.192502)
 - [N3LO-Idaho](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.68.041001)
 - [Malfliet-Tjon](https://www.sciencedirect.com/science/article/pii/0375947469907751)

## Compiling
Download Tic-tac by running
```
git clone https://github.com/seanbsm/Tic-tac.git
```
A makefile example exists in the repository, but there exists currently no way to automatically compile Tic-tac on arbitrary platforms. It is up to the user to compile and link Tic-tac correctly.

## Running the code
There exists a help-functionality built into Tic-tac. Assuming the Tic-tac executable is named `run`, one can run
```
./run --h
```
or 
```
./run -help
```
which displays all available input-arguments to Tic-tac and their function. One can change arguments directly as input on the command-line or write arguments in a `.txt` input file which Tic-tac reads and interprets. If the input-arguments is not of the expected type, for example a float instead of an integer, Tic-tac will stop running.

An example run could be
```
./run Input/input_example.txt Np_WP=32
```
Tic-tac will read `input_example.txt` and change the default input argument values. Non-specified arguments will keep their default value. Tic-tac can read the same input twice, in which case the last read argument will be used during execution. In this example, `Np_WP=32` will rewrite the default value `Np_WP=30` as well as the specified input in `input_example.txt` since it appears last in the command-line.

## How to cite
We prefer that Tic-tac is referred to by name. For reference, use the [original publication](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.106.024001).
