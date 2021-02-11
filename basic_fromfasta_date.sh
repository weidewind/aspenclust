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
for protein in n1pandemic
do
	perl test_global.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --verbose --date &
	perl test_single_sites.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --date &
	perl test_groups_batch.pl --input $input --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --verbose --date &

done
wait
