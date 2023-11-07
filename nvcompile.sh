#!/bin/bash

nvfortran -c hex_utils.f90 hex_main.f90
nvfortran hex_utils.o hex_main.o -o hex_main

time ./hex_main
