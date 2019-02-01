#!/bin/bash

for index in {1..5}; do
	a=$
	b=cat
	cd ~/"TestScripts"
	echo "$a$b KPOINTS MgO
	0
	M
	$index $index $index
	0 0 0" >> KPOINTS$index

	done

