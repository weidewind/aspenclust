#!/bin/bash
while getopts p:t: option
do
                case "${option}" in
                        p) prot=${OPTARG};;
                        t) bigt=${OPTARG};;
                esac
done
number=100
perl test_single_sites.pl --state syn --stattype median --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood  &
perl test_single_sites.pl --state syn --stattype mean --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood  &
perl test_global.pl --state syn --stattype median --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood --verbose  &
perl test_global.pl --state syn --stattype mean --protein $prot --simnumber $number --norm weightnorm --bigtag $bigt --likelihood  --verbose  &
wait
