#!/bin/bash

time ./hex_main 1.90264929971534712693433092186275811e-2 1e-2
mkdir results/1.902649/1e-2
cp *.dat results/1.902694/1e-2/
rm *.dat

time ./hex_main 1.90264929971534712693433092186275811e-2 5e-3
mkdir results/1.902649/5e-3
cp *.dat results/1.902649/5e-3/
rm *.dat

time ./hex_main 1.90264929971534712693433092186275811e-2 1e-3
mkdir results/1.902649/1e-3
cp *.dat results/1.902649/1e-3/
rm *.dat

time ./hex_main 1.90264929971534712693433092186275811e-2 5e-4
mkdir results/1.902649/5e-4
cp *.dat results/1.902649/5e-4/
rm *.dat

time ./hex_main 1.90264929971534712693433092186275811e-2 1e-5
mkdir results/1.902649/1e-5
cp *.dat results/1.902649/1e-5/
rm *.dat
