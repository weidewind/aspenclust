#!/usr/bin/bash

simnum=1000
analysis="all"
state="syn"

for prot in h1 n1 h3 n2 h1pandemic n1pandemic
	do
	perl test_splits.pl --protein ${prot} --state $state --analysis $analysis --stattype mean \
	--splits 4split_mega --input little_ksu --simnumber ${simnum} --verbose  --bigtag test${simnum} &
	done
wait
