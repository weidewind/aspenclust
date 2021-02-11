#!/bin/bash
while getopts t:s:i: option
do
                case "${option}" in
                        t) tag=${OPTARG};;
			s) state=${OPTARG};;
			i) iterations=${OPTARG};;
                esac
done
for protein in h1 h3 n1 n2 
do
	perl test_single_sites.pl --state $state --stattype mean --protein $protein --simnumber $iterations --norm weightnorm --bigtag $tag --likelihood &
	perl test_global.pl --state $state --stattype mean --protein $protein --simnumber $iterations --norm weightnorm --bigtag $tag --likelihood  --verbose  &
	perl test_single_sites.pl --state $state --stattype mean --protein $protein --simnumber $iterations --norm countnorm --bigtag $tag --likelihood  &
	perl test_global.pl --state $state --stattype mean --protein $protein --simnumber $iterations --norm countnorm --bigtag $tag --likelihood  --verbose  &
done
wait
