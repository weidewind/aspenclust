#!/usr/bin/bash

simnum=1000
fakenum=100
state="nsyn"

for prot in h1pdm n1pdm
	do
	echo "fdr_single_sites for $prot"
	perl fdr_ancestors.pl --input knownpass --splits megasplits --state $state --stattype mean --protein $prot --simnumber $simnum --fakenumber $fakenum --bigtag fdrtest_ancestors_${fakenum}_${simnum}
	#echo "grep ss fdr for $prot"
	#python grep_ancestors_fdr.py --fdr /export/home/popova/workspace/aspenclust/output/knownpass/megasplits/fdrtest${fakenum}_${simnum}/${state}/mostlikely/weightnorm --input /export/home/popova/workspace/aspenclust/output/knownpass/megasplits/test10000/${state}/mostlikely/weightnorm --state $state
	done
