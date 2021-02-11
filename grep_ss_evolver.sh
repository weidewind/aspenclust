#!/bin/bash

input="output/evolver_h3/evolvertest/nsyn/mostlikely/weightnorm/"
unlink "${input}sites_results_mean"
for i in {1..99}
	do
	grep "^[0-9]" ${input}${i}_nsyn_sites_mean_statistics >>${input}sites_results_mean
	done
