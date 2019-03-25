#!/bin/bash

# Make sure you specify the B atom of interest when you run the script.
echo Enter the B atom of interest: 
read B_atom
for file in *.vasp; do
	halide_atom=I
	./maxMinDistortionAngles.py ${file} $B_atom $halide_atom
done

