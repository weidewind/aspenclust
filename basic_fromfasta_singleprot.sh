#!/bin/bash
while getopts t:s:i:n:p: option
do
                case "${option}" in
                        t) tag=${OPTARG};;
			s) state=${OPTARG};;
			i) input=${OPTARG};;
			n) number=${OPTARG};;
			p) protein=${OPTARG};;
                esac
done

perl test_global.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --verbose &
perl test_single_sites.pl --input $input --state $state --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag   &
perl test_groups_batch.pl --input $input --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag --verbose &
perl test_alleles.pl --input $input --stattype mean --protein $protein --simnumber $number --norm countnorm --bigtag $tag &
wait
