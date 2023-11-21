#!/bin/bash

cd src/
gfortran -c hex_utils.f90 calc_angle.f90 -fopenmp -pg
#gfortran src/hex_utils.o src/hex_main.o -o src/hex_main -fcheck=bounds -Wall -Wextra -fbacktrace -Og -fopenmp -fopenacc -pg
cd ..
gfortran src/hex_utils.o src/calc_angle.o -o calc_angle -Wall -Wextra -O2 -fopenmp -pg

time ./calc_angle
