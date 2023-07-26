# hex-lattice-symmetry
Studies of symmetry in hexagonal lattices

The Fortran code generates 2 hexagonal plane lattices separated by a distance `z`, and rotates the second one by a defined angle.
The Python code is used for visualization.

Compilation instructions:
`gfortran -c hex_utils.f90 hex_main.f90`
`gfortran hex_utils.o hex_main.o -o hex_main -g`
`./hex_main` or `gdb hex_main` or `valgrind -s ./hex_main`
