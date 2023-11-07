#!/bin/bash

gfortran -c hex_utils.f90 hex_main.f90 -fopenmp -fopenacc -pg
#gfortran hex_utils.o hex_main.o -o hex_main -fcheck=bounds -Wall -Wextra -fbacktrace -Og -fopenmp -fopenacc -pg
gfortran hex_utils.o hex_main.o -o hex_main -Wall -Wextra -O2 -fopenmp -fopenacc -pg

time ./hex_main
