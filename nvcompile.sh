#!/bin/bash

nvfortran -c hex_utils.f90 hex_main.f90 -fopenmp -pg
nvfortran hex_utils.o hex_main.o -o hex_main -Wall -Wextra -O2 -fopenmp -pg

time ./hex_main
