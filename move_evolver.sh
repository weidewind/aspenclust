#!/bin/bash
for i in {1..100}
	do
		cp ../evopoisson/data/evolver_h3_${i}/h3.l.r.newick data/evolver_h3/${i}.ancestor.likelihoods.withinternals.newick
		grep "^>STR" -A 1 ../evopoisson/data/evolver_h3_${i}/h3.all.fa >data/evolver_h3/${i}.nointernals.fa		
		grep "^>N" -A 1 ../evopoisson/data/evolver_h3_${i}/h3.all.fa >data/evolver_h3/${i}.ancestor.fa
	done
