#!/bin/bash
while getopts t:s:i:n: option
do
                case "${option}" in
                        t) tag=${OPTARG};;
			s) state=${OPTARG};;
			i) input=${OPTARG};;
			n) number=${OPTARG};;
                esac
done
for protein in h1 n1 h3 n2 h1pandemic n1pandemic 
do
	perl test_global.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --verbose --fromcsv  &
	perl test_single_sites.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --fromcsv &
	perl test_groups_batch.pl --input $input --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --verbose --fromcsv &

done
wait
