#!/bin/bash

gfortran -c hex_utils.f90 hex_main.f90 -fopenmp
gfortran hex_utils.o hex_main.o -o hex_main -O2 -fopenmp
./hex_main
