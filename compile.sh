#!/bin/bash

cd src/
gfortran -c hex_utils.f90 hex_main.f90 -fopenmp -pg
#gfortran hex_utils.o hex_main.o -o hex_main -fcheck=bounds -Wall -Wextra -fbacktrace -Og -fopenmp -pg
cd ..
gfortran src/hex_utils.o src/hex_main.o -o hex_main -Wall -Wextra -O2 -fopenmp -pg

# angles(1) = 1.90264929971534712693433092186275811e-2
# angles(2) = 1.90403102830304951869456657448196344e-2
# angles(3) = 1.90495757117080288768798873954979561e-2
# angles(4) = 1.91235783233762423239743767417441124e-2
# angles(5) = 1.91517084211538477391489120822979632e-2
# angles(6) = 1.91756277374731560763308567863988621e-2
