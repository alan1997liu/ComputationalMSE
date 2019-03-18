#!/bin/bash

# Make sure you specify the B atom of interest when you run the script.
echo Enter the B atom of interest: 
read B_atom
echo Enter halide atom of interest:
read halide_atom
for file in CONTCAR* ; do
	./general_oct_distortion.py ${file} $B_atom $halide_atom
done

