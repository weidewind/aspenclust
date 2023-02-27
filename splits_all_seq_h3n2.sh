#!/usr/bin/bash

simnum=10000
analysis="all"
state="syn"

for prot in h3 n2
	do
	perl test_splits.pl --protein ${prot} --state $state --analysis $analysis --stattype mean \
	--splits megasplits --input knownpass --simnumber ${simnum} --verbose  --bigtag test${simnum}
	done
