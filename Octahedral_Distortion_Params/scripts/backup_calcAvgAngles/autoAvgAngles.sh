#!/bin/bash

# Make sure you specify the B atom of interest when you run the script.
echo Enter the B atom of interest: 
read B_atom
echo Enter the halide atom of interest: 
read halide_atom
for file in *.vasp; do
	./avgDistortionAngles.py ${file} $B_atom $halide_atom
done

