#!/bin/bash

time ./hex_main 1.90264929971534712693433092186275811e-2
mkdir data/simulation1
cp *.dat data/simulation1/

time ./hex_main 1.90403102830304951869456657448196344e-2
mkdir data/simulation2
cp *.dat data/simulation2/

time ./hex_main 1.90495757117080288768798873954979561e-2
mkdir data/simulation3
cp *.dat data/simulation3/

time ./hex_main 1.91235783233762423239743767417441124e-2
mkdir data/simulation4
cp *.dat data/simulation4/

time ./hex_main 1.91517084211538477391489120822979632e-2
mkdir data/simulation5
cp *.dat data/simulation5/

time ./hex_main 1.91756277374731560763308567863988621e-2
mkdir data/simulation6
cp *.dat data/simulation6/

time ./hex_main 1.91666963330783704217126326302819188e-2
mkdir data/simulation0
cp *.dat data/simulation0/

# angles(1) = 1.90264929971534712693433092186275811e-2
# angles(2) = 1.90403102830304951869456657448196344e-2
# angles(3) = 1.90495757117080288768798873954979561e-2
# angles(4) = 1.91235783233762423239743767417441124e-2
# angles(5) = 1.91517084211538477391489120822979632e-2
# angles(6) = 1.91756277374731560763308567863988621e-2
# angles(7) = 1.91666963330783704217126326302819188e-2
