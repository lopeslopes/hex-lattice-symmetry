# hex-lattice-symmetry
Studies of symmetry in hexagonal lattices

The Fortran code generates 2 hexagonal plane lattices separated by a distance `z`, and rotates the second one by a defined magic angle.
The Python code is used for visualization.

Compilation instructions:  
`gfortran -c hex_utils.f90 hex_main.f90`  
`gfortran hex_utils.o hex_main.o -o hex_main -g`  
`./hex_main` or `gdb hex_main` or `valgrind -s ./hex_main`  

When the lattices are generated with AB stacking, all the points that are found to be overlapping are AB points. When the lattices are generated with AA stacking, all overlapping points are AA or BB.
