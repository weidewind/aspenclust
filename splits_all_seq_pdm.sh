#!/usr/bin/bash

simnum=10000
analysis="all"
state="syn"

for prot in h1pdm n1pdm
	do
	perl test_splits.pl --protein ${prot} --state $state --analysis $analysis --stattype mean \
	--splits megasplits --input knownpass --simnumber ${simnum} --verbose  --bigtag test${simnum}
	done
