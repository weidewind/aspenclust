#!/usr/bin/bash

simnum=1000
fakenum=100
state="syn"

for prot in h3 n2
	do
	echo "fdr_single_sites for $prot"
	perl fdr_single_sites.pl --input knownpass --splits megasplits --state $state --stattype mean --protein $prot --simnumber $simnum --fakenumber $fakenum --bigtag fdrtest${fakenum}_${simnum}
	echo "grep ss fdr for $prot"
	python grep_ss_fdr.py --fdr /export/home/popova/workspace/aspenclust/output/knownpass/megasplits/fdrtest${fakenum}_${simnum}/${state}/mostlikely/weightnorm --input /export/home/popova/workspace/aspenclust/output/knownpass/megasplits/test10000/${state}/mostlikely/weightnorm --state $state
	echo "fdr_alleles for $prot"
	perl fdr_alleles.pl --input knownpass --splits megasplits --state $state --stattype mean --protein $prot --simnumber $simnum --fakenumber $fakenum --bigtag fdrtest${fakenum}_${simnum} 
	echo "grep alleles for $prot"
	python grep_alleles_fdr.py --fdr /export/home/popova/workspace/aspenclust/output/knownpass/megasplits/fdrtest${fakenum}_${simnum}/${state}/mostlikely/weightnorm --input /export/home/popova/workspace/aspenclust/output/knownpass/megasplits/test10000/${state}/mostlikely/weightnorm --state $state
	done
