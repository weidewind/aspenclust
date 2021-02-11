#!/bin/bash

input="output/evolver_h3/evolvertest/nsyn/mostlikely/weightnorm/"
unlink "${input}global_results_mean"
for i in {1..99}
	do
	more ${input}${i}_nsyn_global_mean_statistics | sed -En "s/>/$i/p" >>${input}global_results_mean
	done
