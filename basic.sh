#!/bin/bash
while getopts p:t:s:i: option
do
                case "${option}" in
                        p) protein=${OPTARG};;
                        t) tag=${OPTARG};;
			s) state=${OPTARG};;
			i) iterations=${OPTARG};;
                esac
done
perl test_single_sites.pl --state $state --stattype mean --protein $protein --simnumber $iterations --norm weightnorm --bigtag $tag --likelihood  --verbose &
perl test_global.pl --state $state --stattype mean --protein $protein --simnumber $iterations --norm weightnorm --bigtag $tag --likelihood  --verbose  &
perl test_single_sites.pl --state $state --stattype median --protein $protein --simnumber $iterations --norm weightnorm --bigtag $tag --likelihood  --verbose &
perl test_global.pl --state $state --stattype median --protein $protein --simnumber $iterations --norm weightnorm --bigtag $tag --likelihood  --verbose  &
wait
