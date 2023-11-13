#!/bin/bash

gfortran -c hex_utils.f90 calc_angle.f90 -fopenmp -pg
#gfortran hex_utils.o hex_main.o -o hex_main -fcheck=bounds -Wall -Wextra -fbacktrace -Og -fopenmp -fopenacc -pg
gfortran hex_utils.o calc_angle.o -o calc_angle -Wall -Wextra -O2 -fopenmp -pg

time ./calc_angle
