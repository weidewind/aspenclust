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
	perl test_global.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --verbose &
	perl test_single_sites.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag   &
	perl test_groups_batch.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --verbose >output/${input}/${tag}/groupslog &
	perl test_alleles.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag >output/${input}/${tag}/alleleslog &
done
wait
