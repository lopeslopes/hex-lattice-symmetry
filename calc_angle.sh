#!/bin/bash

gfortran -c /src/hex_utils.f90 /src/calc_angle.f90 -fopenmp -pg
#gfortran /src/hex_utils.o /src/hex_main.o -o /src/hex_main -fcheck=bounds -Wall -Wextra -fbacktrace -Og -fopenmp -fopenacc -pg
gfortran /src/hex_utils.o /src/calc_angle.o -o /src/calc_angle -Wall -Wextra -O2 -fopenmp -pg

time ./src/calc_angle
