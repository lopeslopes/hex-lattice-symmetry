#!/bin/bash

cd src/
gfortran -c hex_utils.f90 hex_main.f90 -fopenmp -pg
#gfortran hex_utils.o hex_main.o -o hex_main -fcheck=bounds -Wall -Wextra -fbacktrace -Og -fopenmp -pg
cd ..
gfortran src/hex_utils.o src/hex_main.o -o hex_main -Wall -Wextra -O2 -fopenmp -pg
