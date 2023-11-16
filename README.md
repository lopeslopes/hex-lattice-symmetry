# hex-lattice-symmetry
Studies of symmetry in hexagonal lattices

The Fortran code generates 2 hexagonal plane lattices separated by a distance `z`, and rotates the second one by a defined magic angle.
The Python code is used for visualization.

Compilation instructions:  
`./compile.sh` for the main program.  
`./calc_angle.sh` to obtain new angles.  
`./nvcompile.sh` for Nvidia Fortran compilator (still not working).
`./remove.sh` to clean .dat files from previous runs.
