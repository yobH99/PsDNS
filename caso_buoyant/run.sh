#!/usr/bin/env bash

# go into data folder and remove all files
cd data
rm -rf *

# go into spectra_2D folder and remove all files
cd ../spectra_2D
rm -rf *

# return to original directory
cd ..

# run the MPI job
mpirun --oversubscribe -n 16 python3 buoyant.py

python3 to_vtk
