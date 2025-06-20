#!/bin/bash


pops=(
	"NonPBO"
	"PBO"
)
chroms=(
	"2RL"
	"3RL"
	"X"
)

for chrom in "${chroms[@]}"
do
	for pop in "${pops[@]}"
	do
		python3 Harun_h12_Feb6.py "$pop" "$chrom" 1000 &
	done
	wait
done


