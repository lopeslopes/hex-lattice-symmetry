#!/bin/bash

time ./hex_main 1.90264929971534712693433092186275811e-2 1e-2
mkdir data/simulation1
cp *.dat data/simulation1/

time ./hex_main 1.90264929971534712693433092186275811e-2 5e-3
mkdir data/simulation2
cp *.dat data/simulation2/

time ./hex_main 1.90264929971534712693433092186275811e-2 1e-3
mkdir data/simulation3
cp *.dat data/simulation3/

time ./hex_main 1.90264929971534712693433092186275811e-2 5e-4
mkdir data/simulation4
cp *.dat data/simulation4/

time ./hex_main 1.90264929971534712693433092186275811e-2 1e-5
mkdir data/simulation5
cp *.dat data/simulation5/
