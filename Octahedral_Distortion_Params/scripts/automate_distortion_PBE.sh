#!/bin/bash
echo Enter in B atom of interest: 
read B_atom
for file in PBE*; do
	halide_atom=I
	./general_oct_distortion.py ${file} $B_atom $halide_atom
done

