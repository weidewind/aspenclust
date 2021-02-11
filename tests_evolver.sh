#!/bin/bash
num=1000
for i in {0..10}
	do
		for j in {1..10} 
		do
			t=$((i*10+j))
			sleep 3
			( perl test_global.pl --stattype mean --protein ${t} --simnumber $num --norm countnorm --input evolver_h3 --bigtag evolvertest --verbose &
			perl test_single_sites.pl --stattype mean --protein ${t} --simnumber $num --norm countnorm --input evolver_h3 --bigtag evolvertest &
                	perl test_global.pl --stattype mean --protein ${t} --simnumber $num --norm weightnorm --input evolver_h3 --bigtag evolvertest --verbose &
	                perl test_single_sites.pl --stattype mean --protein ${t} --simnumber $num --norm weightnorm --input evolver_h3 --bigtag evolvertest ) &

		done
		wait
	done
